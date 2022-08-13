%% set paths: trodes file loader etc
restoredefaultpath;
addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\toolboxes\ums2k_02_23_2012'));
addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared\'));
addpath('C:\Trodes_2-2-3_Windows64\Resources\TrodesToMatlab\');

cd('C:\Data\r204_screening_rec16');
fn = 'r204_screening_rec16.rec';

%% load and reference
ttno = 11;
refno = 9;

out = readTrodesFileContinuous(fn, [ttno 1; ttno 2; ttno 3; ttno 4]);
Fs = out.samplingRate;

out_ref = readTrodesFileContinuous(fn, [refno 1]);
for iC = 1:4
    out.channelData(:, iC) = out.channelData(:, iC) - out_ref.channelData;
end
clear out_ref;
raw = out;

%% filter in spike band
bpFilt = designfilt('bandpassfir', 'StopbandFrequency1', 550, 'PassbandFrequency1', 650, 'PassbandFrequency2', 7000, 'StopbandFrequency2', 8000, 'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', 30000);

for iC = 1:4
    out.channelData(:, iC) = filtfilt(bpFilt, out.channelData(:, iC));
end

%% restrict data for testing purposes
%nSamples = 1e6;
%out.channelData = out.channelData(1:nSamples, :);
%out.timestamps = out.timestamps(1:nSamples);

%% use ums2k spike detector
spk = ss_default_params(Fs);
spk.params.window_size = 1.06; % in ms
spk.params.cross_time = 0.24; % in ms
spk.params.thresh = 4; % SDs above mean
spk.params.max_jitter = 0.5;

clear temp_data; temp_data{1} = out.channelData;
spk = ss_detect(temp_data, spk);
spk = ss_align(spk)

spk.unwrapped_times = spk.unwrapped_times + out.timestamps(1);
clear temp_data;

%% some plotting
nSpks = 100;
spk_idx = nearest_idx3(double(spk.unwrapped_times(1:nSpks)), out.timestamps);

% find out which channel is biggest
wv_temp = sq(max(spk.waveforms, [], 2));
[~, wv_idx] = max(wv_temp, [], 2);

figure(1)
h1 = subplot(211);
spk_idx_lin = sub2ind(size(raw.channelData), spk_idx, wv_idx(1:nSpks));
spk_val = raw.channelData(spk_idx_lin);
plot(out.timestamps(1:spk_idx(end)), raw.channelData(1:spk_idx(end), :));
hold on;
plot(spk.unwrapped_times(1:nSpks), spk_val - 10, '.r', 'MarkerSize', 20);

h2 = subplot(212);
spk_val = out.channelData(spk_idx_lin);
plot(out.timestamps(1:spk_idx(end)), out.channelData(1:spk_idx(end), :));
hold on;
plot(spk.unwrapped_times(1:nSpks), spk_val - 10, '.r', 'MarkerSize', 20);
linkaxes([h1 h2], 'x');

for iSpk = 1:length(spk.unwrapped_times)
    figure(2)
    plot(sq(spk.waveforms(iSpk,:,:))); 
    figure(1)
    xlim([spk.unwrapped_times(iSpk)-0.005 spk.unwrapped_times(iSpk)+0.005])
    
    pause; 
end
%%
wv_peak = max(spk.waveforms, [], 2);
wv_valley = min(spk.waveforms, [], 2);
subplot(331); plot(wv_peak(:,1),wv_peak(:,2),'.k'); axis off;
subplot(332); plot(wv_peak(:,1),wv_peak(:,3),'.k'); axis off;
subplot(333); plot(wv_peak(:,1),wv_peak(:,4),'.k'); axis off;
subplot(334); plot(wv_peak(:,2),wv_peak(:,3),'.k'); axis off;
subplot(335); plot(wv_peak(:,2),wv_peak(:,4),'.k'); axis off;
subplot(336); plot(wv_valley(:,1),wv_peak(:,1),'.r'); axis off;
subplot(337); plot(wv_valley(:,2),wv_peak(:,2),'.r'); axis off;
subplot(338); plot(wv_valley(:,3),wv_peak(:,3),'.r'); axis off;
subplot(339); plot(wv_valley(:,4),wv_peak(:,4),'.r'); axis off;

%%
%spk.waveforms = -double(spk.waveforms) * 100; % invert spikes, and conversion factor hack to be able to see waveforms better in MClust
spk.waveforms = double(spk.waveforms) * 100; % non-inverted spikes...
[~, fname, ~] = fileparts(fn);
fn_out = cat(2, fname, '_TT', num2str(ttno));
save(fn_out, 'spk');

% %% load example header -- need this to write later
% cd('C:\Data\M090\M090-2020-10-02');
% [~, ~, ~, ~, ~, Header] = Nlx2MatSpike('M090-2020-10-02-TT02.ntt', [1 1 1 1 1], 1, 1, []);
% 
% %% try export: data types taken from https://neuralynx.com/_software/NeuralynxDataFileFormats.pdf
% Timestamps = uint64(round(spk.unwrapped_times .* 1e9)); % needs to be in microseconds
% ScNumbers = uint32(ones(size(Timestamps))); % Nlx internal, not used in rest of pipeline
% CellNumbers = uint32(zeros(size(Timestamps))); % Nlx internal, not used in rest of pipeline
% Features = uint32(zeros([8 length(Timestamps)])); % Nlx internal, not used in rest of pipeline
% Samples = int16(permute(-spk.waveforms, [2 3 1]));
% %Samples = zeros(size(Samples));
% Mat2NlxSpike('test5.ntt', 0, 1, [], [1 1 1 1 1], Timestamps, ScNumbers, CellNumbers, Features, Samples, Header); % crashes; setting Samples to zeros fixes that
% 
% %% try to recover -- this runs, but output arguments area all zeros...
% [TimestampsR, ScNumbersR, CellNumbersR, FeaturesR, SamplesR, HeaderR] = Nlx2MatSpike('test5.ntt', [1 1 1 1 1], 1, 1, []);
% 
% %% export simple test data: data types taken from https://neuralynx.com/_software/NeuralynxDataFileFormats.pdf
% Timestamps = uint64(1e9 * [1 2 3]); % needs to be in microseconds
% ScNumbers = uint32([1 1 1]); % Nlx internal, not used in rest of pipeline
% CellNumbers = uint32([0 0 0]); % Nlx internal, not used in rest of pipeline
% Features = uint32(zeros(8, length(Timestamps))); % Nlx internal, not used in rest of pipeline
% Samples = int16(zeros(32, 4, length(Timestamps))); % these are the actual waveforms
% %Header = [];
% Mat2NlxSpike('test2.ntt', 0, 1, [], [1 1 1 1 1], Timestamps, ScNumbers, CellNumbers, Features, Samples, Header);
% 
% %% try to recover -- runs, but recovered timestamps are all zero
% [TimestampsR, ScNumbersR, CellNumbersR, FeaturesR, SamplesR, HeaderR] = Nlx2MatSpike('test2.ntt', [1 1 1 1 1], 1, 1, []);

