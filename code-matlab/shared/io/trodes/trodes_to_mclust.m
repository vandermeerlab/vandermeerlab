%% set paths: trodes file loader etc
restoredefaultpath;
addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\toolboxes\ums2k_02_23_2012'));
%addpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\shared\io\neuralynx');
addpath('C:\Trodes_2-2-3_Windows64\Resources\TrodesToMatlab\');

%%
cd('C:\Data\R204\r204_screening_sq2');
fn = 'r204_screening_sq2.rec';
ttno = 3;

out = readTrodesFileContinuous(fn, [ttno 1; ttno 2; ttno 3; ttno 4]);
Fs = out.samplingRate;

%% 
Wp = [ 700 8000] * 2 / Fs; % pass band for filtering
Ws = [ 500 10000] * 2 / Fs; % transition zone
[N, Wn] = buttord(Wp, Ws, 3, 20); % determine filter parameters
[B, A] = butter(N, Wn); % builds filter

bpFilt = designfilt('bandpassfir', 'FilterOrder', 1000, 'CutoffFrequency1', 600, 'CutoffFrequency2', 9000, 'SampleRate', Fs);

for iC = 1:4
    %out.channelData(:, iC) = filtfilt(bpFilt, out.channelData(:, iC));
    out.channelData(:, iC) = filtfilt(B, A, out.channelData(:, iC));
end

%% restrict data for testing purposes
%nSamples = 1e6;
%out.channelData = out.channelData(1:nSamples, :);
%out.timestamps = out.timestamps(1:nSamples);

%% use ums2k spike detector
spk = ss_default_params(Fs);
spk.params.window_size = 1.06; % in ms
spk.params.cross_time = 0.3; % in ms
spk.params.thresh = 4; % SDs above mean
spk.params.max_jitter = 0.5;

clear temp_data; temp_data{1} = out.channelData;
spk = ss_detect(temp_data, spk);
spk = ss_align(spk)

%%
spk.waveforms = -spk.waveforms * 100; % invert spikes, and conversion factor hack to be able to see waveforms better in MClust
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

