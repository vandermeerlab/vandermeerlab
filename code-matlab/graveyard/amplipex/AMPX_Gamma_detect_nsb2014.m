function [low_gamma_iv, high_gamma_iv, random_lg_iv, random_hg_iv, noise_iv, channel_of_interest] = AMPX_Gamma_detect_nsb2014 (data, varargin)
% AMPX_Gamma_detect_nsb
%
% This function will use the new nsb2014 codebase to extract the high and
% low gamma events outside of the periods with a high level of noise in two
% bands
%
% INPUTS
% - data [struct] contains, all the channels as well as a header and tvec.  The data should have already been decimated, DC and artifacts already removed
%
% optional:
% - low_gamma_band [1x2 array] as "[low_pass high_pass]" the detection range for low gamma events
% - high_gamma_band [1x2 array] as "[low_pass high_pass]" the detection range for high gamma events
% - noise_band [1x2 array] as "[low_pass high_pass]" the detection range for lower noise artifacts
% - noise_band [1x2 array] as "[low_pass high_pass]" the detection range for higher frequency noise artifacts
%
% OUTPUTS:
% low_gamma_iv [struct] - contains:
%     - tstart: event start time
%     - tend: event end time
%     -

%% input Variables

if isfield(data, 'dc_remove')==0 || isfield(data, 'artifacts_removed')==0
    warning('data has not been preprocessed properly.  Either the DC or the artifacts have not been removed')
    fprintf('\nAre you sure you want to proceed? \nPress any key to proceed...\n')
else
    fprintf('Data has been preprocessed and is ready for event detection\n')
end

if exist('data', 'var') == 0
    load post_data_art_rem_DC_dec10.mat
end
low_gamma_band = [40 55];  % default 40-55
high_gamma_band = [70 85]; % default 70-85
noise_band = [110 225];   % should be in the chewing/scratching range (110-225)
noise_band2 = [225 400];  % should be high ~225-400Hz
channel_of_interest = 1; % which channel is used for detection
display = 'off';           % plots the events and lets your see each event
extract_varargin
%     cfg.fc = {ExpKeys.goodGamma_vStr{1}};
%     data = ft_read_neuralynx_interp(cfg.fc);
%     data.label = {'vStr','mPFC'};
%     data.hdr.Fs = data.fsample; % should be set by ft_read_neuralynx_interp
%% ensure that there are no extreme artifacts.
data_std = std(data.channels{channel_of_interest});
data.channels{channel_of_interest}(data.channels{channel_of_interest} > 3*data_std) =0;
data.channels{channel_of_interest}(data.channels{channel_of_interest} < -3*data_std) =0;
% for ii = 1:length(data.channels{62})
% if data.channels{channel_of_interest}(ii) > 5*data_std || data.channels{channel_of_interest}(ii) < -5* data_std;
% data.channels{62}(ii) = 0;
% end
% end
%% prepare data for event detection
cfg = [];
if length(data.tvec) ==1
    data.tvec = data.tvec'; % transpose the data (not the tvec) or else you will get an erroe in the
    disp('data.tvec needs to be transposed')
end
data.channels{channel_of_interest} = data.channels{channel_of_interest}';
data.cgf = cfg;
data_tsd = tsd(data.tvec,data.channels{channel_of_interest});
data_tsd.cfg.hdr{1}.SamplingFrequency = data.hdr.Fs;
data_tsd.data(isnan(data_tsd.data)) = 0;
data_tsd.label = num2str(channel_of_interest);
cfg = []; cfg.decimateFactor = 1; % reduce sampling frequency to speed things up
csc = decimate_tsd(cfg,data_tsd);
%% chewing artifact detection
cfg = []; cfg.f = noise_band;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
cscF = zscore_tsd(cscF); % normalize

cfg = [];
cfg.method = 'raw';
cfg.threshold = 3.5;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.minlen = 0.01; % minimum interval length

noise_iv = TSDtoIV(cfg,cscF); % detect intervals where z-scored power is > 3 SD above mean

cfg = []; cfg.f = noise_band2;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
cscF = zscore_tsd(cscF); % normalize

cfg = [];
cfg.method = 'raw';
cfg.threshold = 3.5;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.minlen = 0.01; % minimum interval length

noise_iv = UnionIV([],noise_iv,TSDtoIV(cfg,cscF));

cfg = []; cfg.d = [-0.1 0.1];
noise_iv = ResizeIV(cfg,noise_iv);
%%
if strcmp(display, 'on')
    PlotTSDfromIV([],noise_iv,csc);
end
%% low gamma detection
cfg = []; cfg.f = low_gamma_band;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
cscF = zscore_tsd(cscF); % normalize

cfg = [];
cfg.method = 'raw';
cfg.threshold = 2;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.025; % merge events closer than this
cfg.minlen = 0.06; % minimum interval length

lg_iv = TSDtoIV(cfg,cscF); % detect intervals where z-scored power is > 3 SD above mean
lg_iv = DifferenceIV([],lg_iv,noise_iv);
fprintf(['Events detected: ' num2str(length(lg_iv.tstart)) '\n'])
%% low gamma low thres detection
cfg = []; cfg.f = low_gamma_band;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
cscF = zscore_tsd(cscF); % normalize

cfg = [];
cfg.method = 'raw';
cfg.threshold = .75;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.025; % merge events closer than this
cfg.minlen = 0.06; % minimum interval length

lg_low_iv = TSDtoIV(cfg,cscF); % detect intervals where z-scored power is > 3 SD above mean
lg_low_iv = DifferenceIV([],lg_low_iv,noise_iv);
fprintf(['Events detected: ' num2str(length(lg_low_iv.tstart)) '\n'])
%% random time points
r = sort(randi([1 length(data.tvec)], 1, length(lg_iv.tstart)));
for ii = length(r):-1:1
    ran_lg_iv.tstart(ii,1) = data.tvec(r(ii));
    ran_lg_iv.tend(ii,1) = ran_lg_iv.tstart(ii) + lg_iv.tend(ii)-lg_iv.tstart(ii);
end
ran_lg_iv.cfg = lg_iv.cfg;
%%
if strcmp(display, 'on')
    close all
    PlotTSDfromIV([],lg_iv,csc);
    for iEvent = 1:length(lg_iv.tstart)
        set(gca,'XLim',[lg_iv.tstart(iEvent)-0.5 lg_iv.tstart(iEvent)+.5]);
        text(lg_iv.tstart(iEvent), 100, num2str(iEvent))
        pause
    end
end
%% high gamma detection
cfg = []; cfg.f = high_gamma_band;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
cscF = zscore_tsd(cscF); % normalize

cfg = [];
cfg.method = 'raw';
cfg.threshold = 2;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.025; % merge events closer than this
cfg.minlen = 0.06; % minimum interval length

hg_iv = TSDtoIV(cfg,cscF); % detect intervals where z-scored power is > 3 SD above mean
hg_iv = DifferenceIV([],hg_iv,noise_iv);
fprintf(['Threshold: ' num2str(cfg.threshold) '   ' 'Min Length:  ' num2str(cfg.minlen) '\n' 'Events detected: ' num2str(length(hg_iv.tstart)) '\n'])
%% high gamma lower detection
cfg = []; cfg.f = high_gamma_band;
cscF = FilterLFP(cfg,csc); % filter
cscF = LFPpower([],cscF); % obtain instantaneous signal power (Hilbert transform)
cscF = zscore_tsd(cscF); % normalize

cfg = [];
cfg.method = 'raw';
cfg.threshold = .75;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.025; % merge events closer than this
cfg.minlen = 0.06; % minimum interval length

hg_low_iv = TSDtoIV(cfg,cscF); % detect intervals where z-scored power is > 3 SD above mean
hg_low_iv = DifferenceIV([],hg_low_iv,noise_iv);
fprintf(['Threshold: ' num2str(cfg.threshold) '   ' 'Min Length:  ' num2str(cfg.minlen) '\n' 'Events detected: ' num2str(length(hg_low_iv.tstart)) '\n'])

%% plot the events
if strcmp(display, 'on')
    close all
    PlotTSDfromIV([],hg_iv,csc);
    for iEvent = 1:length(hg_iv.tstart)
        set(gca,'XLim',[hg_iv.tstart(iEvent)-0.5 hg_iv.tstart(iEvent)+.5]);
        pause
    end
    close all
end
%%
r = sort(randi([1 length(data.tvec)], 1, length(hg_iv.tstart)));
for ii = length(r):-1:1
    ran_hg_iv.tstart(ii,1) = data.tvec(r(ii));
    ran_hg_iv.tend(ii,1) = ran_hg_iv.tstart(ii) + hg_iv.tend(ii)-hg_iv.tstart(ii);
end
ran_hg_iv.cfg = hg_iv.cfg;

%% compile the random events such that they do not over lap with high or low gamma
ran_iv.tstart = [ran_hg_iv.tstart; ran_lg_iv.tstart];
ran_iv.tend = [ran_hg_iv.tend; ran_lg_iv.tend];
ran_iv.cfg = lg_iv.cfg;
fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before low Gamma DifferenceIV \n'])

pre_diff = length(ran_iv.tend);
ran_iv = DifferenceIV([], ran_iv, lg_low_iv);
disp([num2str(length(ran_iv.tstart)) '# of ran events after low Gamma DifferenceIV ' ])
post_diff = length(ran_iv.tend);
fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])

fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before high Gamma DifferenceIV \n'])
pre_diff = length(ran_iv.tend);
ran_iv = DifferenceIV([], ran_iv,hg_low_iv);
disp([num2str(length(ran_iv.tstart)) '# of ran events after high Gamma DifferenceIV ' ])
post_diff = length(ran_iv.tend);
fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])
%% loop this until there is only a group of ran_ivs that are the same length as the lg and Hg but not overlapping with either
loop_num = 0;
checks = 10000:5000:400000;
rng shuffle
while length(ran_iv.tstart) ~= (length(lg_iv.tstart)+length(hg_iv.tstart))
    if ismember(loop_num, checks)
        hg_iv.tstart = hg_iv.tstart(1:end-5);
        hg_iv.tend = hg_iv.tend(1:end-5);
        lg_iv.tstart = lg_iv.tstart(1:end-5);
        lg_iv.tend = lg_iv.tend(1:end-5);
        rng shuffle
        disp(['Ran not hitting at ' num2str(loop_num)])
    end
    loop_num = 1+loop_num;
    %     fprintf(['\n\nStarting loop #' num2str(loop_num) '\n\n'])
    
    r = randi([1 length(data.tvec)], 1, ((length(lg_iv.tstart) + length(hg_iv.tstart))- length(ran_iv.tstart)));
    ran_iv_app = [];
    for ii = length(r):-1:1
        ran_iv_app.tstart(ii,1) = data.tvec(r(ii));
        ran_iv_app.tend(ii,1) = ran_iv_app.tstart(ii) + hg_low_iv.tend(ii)-hg_low_iv.tstart(ii);
    end
    ran_iv_app.cfg = ran_iv.cfg;
    pre_diff = length(ran_iv.tend);
    ran_iv = DifferenceIV([], ran_iv, ran_iv_app);
    post_diff = length(ran_iv.tend);
    %     fprintf(['\n' num2str(pre_diff-post_diff) '  Events removed ran_iv VS ran_iv_app\n'])
    if (pre_diff-post_diff) ==0
        ran_iv.tstart = [ran_iv.tstart; ran_iv_app.tstart];
        ran_iv.tend = [ran_iv.tend; ran_iv_app.tend];
    else
        continue
    end
    %     fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before low Gamma DifferenceIV \n'])
    pre_diff = length(ran_iv.tend);
    ran_iv = DifferenceIV([], ran_iv, lg_low_iv);
    %     disp([num2str(length(ran_iv.tstart)) '# of ran events after low Gamma DifferenceIV ' ])
    post_diff = length(ran_iv.tend);
    %     fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])
    
    %     fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before high Gamma DifferenceIV \n'])
    pre_diff = length(ran_iv.tend);
    ran_iv = DifferenceIV([], ran_iv,hg_low_iv);
    %     disp([num2str(length(ran_iv.tstart)) '# of ran events after high Gamma DifferenceIV ' ])
    post_diff = length(ran_iv.tend);
    %     fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])
    %     disp(['Loop Number ' num2str(loop_num)])
end

fprintf(['\n Finished at loop #' num2str(loop_num) '\n\n'])

%% split the random values into two arrays for low and high gamma
ran_lg_iv = []; ran_hg_iv = [];
r_lg = randperm(length(lg_iv.tstart));
ran_lg_iv.tstart = ran_iv.tstart(r_lg);
ran_lg_iv.tend = ran_iv.tend(r_lg);

ran_hg_iv.tstart = ran_iv.tstart((ismember(ran_iv.tstart, ran_lg_iv.tstart))==0);
ran_hg_iv.tend= ran_iv.tend((ismember(ran_iv.tend, ran_lg_iv.tend))==0);

ran_lg_iv.tstart = sort(ran_lg_iv.tstart);
ran_lg_iv.tend = sort(ran_lg_iv.tend);
ran_hg_iv.tstart = sort(ran_hg_iv.tstart);
ran_hg_iv.tend = sort(ran_hg_iv.tend);

%% plot the events
if strcmp(display, 'on')
    close all
    PlotTSDfromIV([],ran_lg_iv,csc);
    for iEvent = 1:length(ran_lg_iv.tstart)
        set(gca,'XLim',[ran_lg_iv.tstart(iEvent)-0.5 ran_lg_iv.tstart(iEvent)+.5]);
        pause
    end
    close all
end
%% output the events.
high_gamma_iv = hg_iv;
low_gamma_iv = lg_iv;
random_lg_iv = ran_lg_iv;
random_hg_iv = ran_hg_iv;

end
