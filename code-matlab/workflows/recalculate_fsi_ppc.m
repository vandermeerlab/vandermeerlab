%% Boiler plate to recalculate spectral measures for any FSI or MSN
% setup

clear;
%load significant cells
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
load('./significance.mat');

% Change the below line to msn or fsi to choose cell accordingly 
label = sig_fsi{randi(length(sig_fsi),1)};

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
%% Calculate unsampled values
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
% Sanity check to ensure that redefine trial has the same number of
% spikes as restrict()
this_flag = false;
spk_count1 = length(restrict(sd.S, rt_iv).t{iC});
spk_count2 = 0;
for iT = 1:length(near_data.trial)
   spk_count2 = spk_count2 + sum(near_data.trial{iT}(2,:)); 
end
od.spk_count = spk_count2;
od.flag_tooFewSpikes = false;
if spk_count1 ~= spk_count2
    this_flag = true;
    warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in near-Reward trials', ...
        spk_count1, spk_count2)
end

% Calculate and save STA
this_sta = ft_spiketriggeredaverage(cfg_ft, near_data);
od.sta_time = this_sta.time;
od.unsampled_sta = this_sta.avg(:,:)';
od.flag_unequalSpikes = this_flag;

% Calculate and save STS
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
od.freqs = this_sts.freq;
od.unsampled_sts = nanmean(sq(abs(this_sts.fourierspctrm{1})));
od.flag_nansts = this_flag;

% Calculate and save PPC (use both ppc0 and ppc2)
cfg_ppc               = [];
cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
cfg_ppc.spikechannel  = this_sts.label;
cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
cfg_ppc.avgoverchan   = 'weighted';
cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
% Display warning to show that there were Nans in this calculation
if ~isempty(find(isnan(this_ppc.ppc0),1))
    warning('Cell %s has nans in the unsampled ppc0',sd.S.label{iC});
end
od.unsampled_ppc0 = this_ppc.ppc0;

cfg_ppc.method = 'ppc2';
this_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
% Display warning to show that there were Nans in this calculation
if ~isempty(find(isnan(this_ppc.ppc2),1))
    warning('Cell %s has nans in the unsampled ppc2',sd.S.label{iC});
end
od.unsampled_ppc2 = this_ppc.ppc2;

%% Set parameters for subsampling
s_factors = [2, 5, 10, 25, 50]; % factor by which spike count is divided in a given round of subsampling
n_subs = 250; % Rounds of subsampling
all_scount = od.spk_count;
od.s_factors = s_factors;
od.n_subs = n_subs;
% Random subsampling
od.rsample = cell(1,length(s_factors));
for iS = 1:length(s_factors)
    this_scount = round(all_scount/s_factors(iS));
    temp_sts = zeros(n_subs, length(od.unsampled_sts));
    temp_ppc = zeros(n_subs, length(od.unsampled_ppc0));
    for iR = 1:n_subs
        sub_idx = randsample(1:all_scount, this_scount);
        sub_idx = sort(sub_idx); sub_sts = this_sts;
        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
        temp_sts(iR,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
        cfg_ppc.method = 'ppc0';
        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
        temp_ppc(iR,:) = sub_ppc.ppc0';
    end
    od.rsample{iS}.sts = mean(temp_sts,1);
    od.rsample{iS}.ppc0 = mean(temp_ppc,1);
end

% Contiguous subsampling
od.csample = cell(1,length(s_factors));
for iS = 1:length(s_factors)
    this_scount = round(all_scount/s_factors(iS));
    temp_sts = zeros(n_subs, length(od.unsampled_sts));
    temp_ppc = zeros(n_subs, length(od.unsampled_ppc0));
    for iR = 1:n_subs
        sub_idx = randi(all_scount, 1);
        sub_idx = [sub_idx:sub_idx+this_scount-1];
        sub_idx(sub_idx > all_scount) = sub_idx(sub_idx > all_scount) - all_scount; %circular shifting indices that exceed total spike count
        sub_idx = sort(sub_idx); sub_sts = this_sts;
        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
        temp_sts(iR,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
        cfg_ppc.method = 'ppc0';
        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
        temp_ppc(iR,:) = sub_ppc.ppc0';
    end
    od.csample{iS}.sts = mean(temp_sts,1);
    od.csample{iS}.ppc0 = mean(temp_ppc,1);
end
%% Save recalculated values
save(strcat('D:\RandomVstrAnalysis\temp\',label, '_ppc_test.mat'), 'od')
%% Plotting
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
fig = figure('WindowState', 'maximized');
for iF = 1:length(f_list)
    ax = subplot(1, length(f_list), iF);
    hold on
    f_idx = find(round(od.freqs) >= f_list{iF}(1) & ...
        round(od.freqs) <= f_list{iF}(2));
    u0 = mean(od.unsampled_ppc0(f_idx));
    u2 = mean(od.unsampled_ppc2(f_idx));
    % take max PPC in the range
%     plot(od.s_factors, u0);
%     plot(od.s_factors, u2);
    plot(od.s_factors, abs(cellfun(@(x) mean(x.ppc0(f_idx)), od.rsample) - u0));
    plot(od.s_factors, abs(cellfun(@(x) mean(x.ppc0(f_idx)), od.csample) - u0));
    ax.YAxis.Exponent = 0;
    xticks(od.s_factors)
    xlabel('Subsampling Factor');
    legend({'Random', 'Contiguous'}, 'FontSize', 10, 'Location','northwest');
    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
end
sgtitle(sprintf('FSI, Total Spikes: %d, Rounds of Subsampling: %d', od.spk_count, od.n_subs))
WriteFig(fig,strcat('D:\RandomVstrAnalysis\temp\',label, '_PPCTest'),1);
close;

