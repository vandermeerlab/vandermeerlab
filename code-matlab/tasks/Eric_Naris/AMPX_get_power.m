function power = AMPX_get_power(cfg_in, evts, data_ft)
%% AMPX_get_power: extracts the power within a specified frequency range
%
%
%   INPUTS:
%      - cfg_in [struct]: contains parameters for the ft_freqanalysis
%      function as well as the frequency range of interest.  
%      - evts [struct]: contains all the ivs for low, high gamma events and
%      thier matched control ivs
%      - data_ft [struct]: FieldTrip data format that has not bee
%      trilified. 
%




%% parameters
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