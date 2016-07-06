function power = AMPX_get_power(cfg_in, data_ft, ExpKeys)
%% AMPX_get_power: extracts the power within a specified frequency range
%
%
%   INPUTS:
%      - cfg_in [struct]: contains parameters for the ft_freqanalysis
%      function as well as the frequency range of interest.
%      - data_ft [struct]: FieldTrip data format that has not bee
%      trilified.
%      - ExpKeys [struct]: contains all the experimental parameters like
%      broken channels, the probe type...
%    OUTPUTS:
%      - power [struct]: contains the power within the band of interest
%        across all channels for each event in evts. Contains:
%               - power [8 x 8 x nEvents]: mean power across event in FOI
%               - cfg [struct]: contains all the cfg values used
%      - power_avg [8x8]: mean or median (cfg.avg_method) power across all
%        events
%
%
%% parameters
cfg_def              = [];
cfg_def.output       = 'pow';
cfg_def.channel      = 'all';
cfg_def.method       = 'mtmfft';
cfg_def.taper        = 'hanning';
cfg_def.foi          = 1:1:500;
cfg_def.pad          = 2;
cfg_def.trials       = 1:length(data_ft.trial); %[evt.post.good_high_gamma];
cfg_def.keeptrials   = 'yes'; % need this for stats later
cfg_def.t_ftimwin    = 5./cfg_def.foi;  % 20 cycles per time window
cfg_def.toi          = 0:0.05:1; % pre-nosepoke baseline
cfg_def.avg_method   = 'mean'; % can be median uses "nanmean" or "nanmedian"
% cfg_def.ft_version = ft_version;
 
cfg = ProcessConfig(cfg_def, cfg_in);
%% compute the power from the ft_freqanalysis function using the cfg parameters
 
TFR = ft_freqanalysis(cfg, data_ft);
 
%% get the power for each event and extract power in the freqency of interest
power_distrib = cell(length(cfg.trials),1);
for itrial = 1:length(cfg.trials)
    for ii =length(TFR.label):-1:1
        if strcmp(cfg.avg_method, 'mean')
            power_distrib{itrial}(1, ii) = nanmean(TFR.powspctrm(itrial,ii, nearest(TFR.freq,cfg.freq(1)):nearest(TFR.freq,cfg.freq(2))));  % finds band of interest
        elseif strcmp(cfg.avg_method, 'median')
            power_distrib{itrial}(1, ii) = nanmedian(TFR.powspctrm(itrial,ii, nearest(TFR.freq,cfg.freq(1)):nearest(TFR.freq,cfg.freq(2))));  % finds band of interest
            
        end
    end
    % nan the bad channels
    for bad = ExpKeys.BadChannels;
        for ii = 1:64
            if ismember(ExpKeys.Probe_layout(ii), bad)
                power_distrib{itrial}(1,ii) = nan;
            end
        end
    end
    power.power_distrib{itrial} = reshape(power_distrib{itrial}, 8,8);  % convert to 8x8 layout.
end
 
%% get the average power across all events in in the 8x8 layout
all_trials_3d = cat(3,power.power_distrib{1,:});
power.power_distrib_avg = mean(all_trials_3d,3);

%% keep the config files with the parameters
power.cfg  = cfg;
power.cfg.date = datestr(date, 'yyyy-mm-dd');