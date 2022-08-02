function [power] = AMPX_pow_distrib (trial_data, varargin)
% function [power_distrib] = AMPX_pow_distrib (trial_data, varargin)
% Takes in trialified data and extracts the power within a frequency range
% of interest as specified by min_freq/max_freq using the Fieldtrip
% spectrogram data obtained using ft_freqanalysis.  This is done on a trial
% by trial basis
%
%To Do Prior to this function:
%   ensure that the ExpKeys is complete and contains the channel layout and
%   the bad channels.  The Expkeys script will be called and used
%   throughout this function.  
% 
%INPUTS: 
%   trial_data: trialified data in the fieldtrip format (trial_data).  This will include
%   a header, tvec, trial data (nchan x length)
%   trials_def 
%
% OUTPUTS: 
%   power is a structure containing:
%   -power_distrib = [1 x nTrials] cell array with each element being the
%    power for each channel
%   -power_distrib_avg = [8 x 8] matrix representing the channels in the
%    probe layout 
%   -params: this contains all the informaiton from the AMPX_pow_distrib
%    variables such as the frequencies used.
%   

% EC v1.0 12/09/2014
%% Define the variables
run(FindFile('*Keys.m'))

min_freq = 45; % low gamma [45 65], high gamma [70 85], high_freq [450 490]
max_freq = 65;

extract_varargin;

%% specgram 
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:500;
cfg.pad          = 2;
cfg.trials       = 1:length(trial_data.trial); %[evt.post.good_high_gamma];
cfg.keeptrials   = 'yes'; % need this for stats later
cfg.t_ftimwin    = 5./cfg.foi;  % 20 cycles per time window
cfg.toi          = 0:0.05:1; % pre-nosepoke baseline
TFR = ft_freqanalysis(cfg, trial_data);

%% extract the power spectrums for each channel across each event
power_distrib = cell(length(cfg.trials),1); 
for itrial = 1:length(cfg.trials)
    for ii =length(TFR.label):-1:1
        power_distrib{itrial}(1, ii) = nanmean(TFR.powspctrm(itrial,ExpKeys.Probe_layout(ii),nearest(TFR.freq,min_freq):nearest(TFR.freq,max_freq)));  % finds the high Gamma power
    end
    % nan the bad channels
    for bad = ExpKeys.BadChannels;
        for ii = 1:64
            if ExpKeys.Probe_layout(ii) == bad
                power_distrib{itrial}(1,ii) = nan;
            end
        end
    end
    power.power_distrib{itrial} = reshape(power_distrib{itrial}, 8,8); 
end


%% Average across all trials.
all_trials_3d = cat(3,power.power_distrib{1,:});
power.power_distrib_avg = mean(all_trials_3d,3);
%% add the paramters to the power structure
power.params.freq = [min_freq max_freq];
power.params.window_size = ((trial_data.cfg.trl(1,2)/trial_data.fsample)-(trial_data.cfg.trl(1,1)/trial_data.fsample))/2;


