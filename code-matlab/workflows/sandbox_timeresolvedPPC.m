%% Boiler plate to recalculate spectral measures for any FSI or MSN

%% setup

clear;
%load significant cells
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
load('./significance.mat');

% Change the below line to msn or fsi to choose cell accordingly 
label = sig_fsi{randi(length(sig_fsi),1)};
label = sig_fsi{2};

% generate correct path on the basis of the label
toks = strsplit(label, '-');
path = strcat('E:\ADRLabData\', toks{1}, '\', strjoin(toks(1:4),'-'));
cd(path);

% Load CSC, ExpKeys and the cell
LoadExpKeys;

if isfield(ExpKeys,'goodGamma_vStr')
    cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
elseif isfield(ExpKeys, 'goodGamma')
    cfg = []; cfg.fc = ExpKeys.goodGamma;
else
    error('Couldn''t find LFP field name.');
end

csc = LoadCSC(cfg); csc.data = csc.data-nanmean(csc.data); % could locdetrend to improve STA estimate
% Note that Fieldtrip will interpoalate data to avoid gaps. Include
% sanity tests to ensure STA/STS segments don't include these gaps.
ft_csc = ft_read_neuralynx_interp(cfg.fc);

% hacky code to separate out only ontrack-data
temp_tvec = [0:length(ft_csc.time{1})-1];
temp_offset = (double(ft_csc.hdr.LastTimeStamp)/1e6 - double(ft_csc.hdr.FirstTimeStamp)/1e6)/(length(temp_tvec) - 1);
temp_tvec = temp_tvec * temp_offset;
 
% Modify ft_csc
ft_csc.time{1} = temp_tvec;
ft_csc.fsample = 1/temp_offset;

% Ensure that the new data doesn't have any Nans in it
temp_tvec = temp_tvec + double(ft_csc.hdr.FirstTimeStamp)/1e6;
temp_start = nearest_idx3(ExpKeys.TimeOnTrack, temp_tvec);
temp_end = nearest_idx3(ExpKeys.TimeOffTrack, temp_tvec);
cfg_onTrack.begsample = temp_start;
cfg_onTrack.endsample = temp_end;
data_onTrack = ft_redefinetrial(cfg_onTrack, ft_csc);

if ~isempty(find(isnan(data_onTrack.time{1}),1))
    warning('On Track data for %s has gaps', cfg.fc{1});
end
%%
iC = 1; iM = 1; % Hard-coding this to make copy-pasting easier

cfg = []; cfg.fc = {strcat(label, '.t')};
sd.S = LoadSpikes(cfg);
% Restrict spikes to only OnTrack
sd.S = restrict(sd.S, iv(ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack));
% Read the spike files into field trip format
sd.S.ft_spikes = ft_read_spike(sd.S.label{1});

% Calculate and save STA
cfg_ft.timwin = [-0.5 0.5];
cfg_ft.spikechannel = sd.S.ft_spikes(iC).label{1};
cfg_ft.channel = ft_csc.label(1);
this_data = ft_appendspike([], ft_csc, sd.S.ft_spikes(iC));

% Block of code to divide recordings session into trials

% restrict spikes to a timeWindow of +/-5 seconds around the reward  
rt1 = getRewardTimes();
rt1 = rt1(rt1 > ExpKeys.TimeOnTrack);
rt2 = getRewardTimes2();
rt2 = rt2(rt2 > ExpKeys.TimeOnTrack);
% Sometimes (in R117-2007-06-12, for instance) getRewardTimes() returns
% times that are spaced out less than 5 sec apart (possibly erroneus). 
% Getting rid of such reward times to maintain consistency
rt_dif = diff(rt1);
rt_dif = find(rt_dif <= 5);
valid_rt1 = true(length(rt1),1);
valid_rt2 = true(length(rt2),1);
for i = 1:length(rt_dif)
    valid_rt1(rt_dif(i)) = false;
    valid_rt1(rt_dif(i)+1) = false;
    valid_rt2(rt2 >= rt1(rt_dif(i)) & rt2 <= rt1(rt_dif(i)+2)) = false;
end
% Sometimes (in R119-2007-07-05, for instance) getRewardTimes2() returns
% times that are spaced out less than 5 sec apart (possibly erroneus). 
% Getting rid of such reward times to maintain consistency
rt_dif = diff(rt2);
rt_dif = find(rt_dif <= 5);
for i = 1:length(rt_dif)
    valid_rt2(rt_dif(i)) = false;
    valid_rt2(rt_dif(i)+1) = false;
    valid_rt1(rt1 >= rt2(rt_dif(i)-1) & rt1 <= rt2(rt_dif(i)+1)) = false;
end
rt1 = rt1(valid_rt1);
rt2 = rt2(valid_rt2);
% Sometimes (in R117-2007-06-17, for instance) the second reward is
% triggered but not the first one in the last trial
if length(rt1) ~= length(rt2)
    rt1 = rt1(1:end-1);
end
% Sanity check to make sure that rt2 is always triggered after rt1
keep = (rt1 <= rt2);
rt1 = rt1(keep);
rt2 = rt2(keep);

% For near reward_trials  
w_start = rt1 - 5;
w_end =  rt2 + 5;
% Last trial time shouldn't exceed Experiment end time
w_end(end) = min(w_end(end), ExpKeys.TimeOffTrack);
% Sorting makes it wonky in some cases (in R119-2007-07-06),so 
% only keep trials that are not outliers
keep = ~isoutlier(w_end - w_start, 'median');
w_start = w_start(keep);
w_end = w_end(keep);
rt_iv = iv(w_start, w_end);

% Break down data into near trials
temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;     
temp_start = nearest_idx3(rt_iv.tstart, temp_tvec);
temp_end = nearest_idx3(rt_iv.tend, temp_tvec);
cfg_near_trials.trl = [temp_start, temp_end, zeros(size(temp_start))];
near_data = ft_redefinetrial(cfg_near_trials, this_data);
spk_count2 = 0;
near_trialwise_spk_count = zeros(1, length(near_data.trial));
for iT = 1:length(near_data.trial)
   near_trialwise_spk_count(iT) = sum(near_data.trial{iT}(2,:));
   spk_count2 = spk_count2 + near_trialwise_spk_count(iT);
end

 % Extract All Spike IDs
trial_wise_spike = cell(1, length(near_data.trial));
last_spk_ct = 0;
for iT = 1:length(near_data.trial)
    this_spk_ct = near_trialwise_spk_count(iT);
    trial_wise_spike{iT} = last_spk_ct + 1 : last_spk_ct + this_spk_ct;
    last_spk_ct = last_spk_ct + this_spk_ct;
end

% Calculate and save STS for all spikes
cfg_sts.method = 'mtmconvol';
cfg_sts.foi = 1:1:100;
cfg_sts.t_ftimwin = 5./cfg_sts.foi;
cfg_sts.taper = 'hanning';
cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
cfg_sts.channel = near_data.label{1};
this_sts = ft_spiketriggeredspectrum(cfg_sts, near_data);
this_flag = false;
% Display warning to show that there were Nans in this calculation
if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
    this_flag = true;
    warning('Cell %s has nans in its STS',sd.S.label{iC});
end
od.fsi_res.near_spec{iM}.freqs = this_sts.freq;
od.fsi_res.near_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
od.fsi_res.near_spec{iM}.flag_nansts = this_flag;

% Calculate and save trialwise_ppc and use subsampled measures
% for near Reward Trials. However, if an FSI has fewer spikes
% than at least one co-recorded MSN, do NOT
% subsample and include all the spikes

% Generate trialwise_ppc
this_tw_ppc = nan(length(near_data.trial), length(od.fsi_res.near_spec{iM}.freqs));
cfg_ppc               = [];
cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
cfg_ppc.spikechannel  = this_sts.label;
cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
cfg_ppc.avgoverchan   = 'weighted';
this_flag = false;

for iT = 1:length(near_data.trial)
    sub_idx = trial_wise_spike{iT};
    sub_sts = this_sts;
    sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
    sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
    sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
    sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
    this_tw_ppc(iT,:) = sub_ppc.ppc0;
end

%% Trying timesolved ppc
cfg_ppc2 = cfg_ppc;
cfg_ppc2.latency = 'maxperiod';
cfg_ppc2.timwin = [0 2];
cfg_ppc2.winstepsize = 1;

for iT = 1:length(near_data.trial)
    sub_idx = trial_wise_spike{iT};
    sub_sts = this_sts;
    sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
    sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
    sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
    sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc2, sub_sts);
    q = squeeze(sub_ppc.nspikes);
    dummy = 1;
end