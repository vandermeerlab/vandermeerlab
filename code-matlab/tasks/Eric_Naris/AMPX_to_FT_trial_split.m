function data_trl = AMPX_to_FT_trial_split(cfg_in,evts_in, data_remap_FT)
%% AMPX_trial_split_data: converts the iv times for the gamma events to
% data trials in the FieldTrip format
%
%   Inputs:
%     - data_remap [struct]: data in the FT format after the channels
%     have been remapped.
%     - evts_in [struct]: contains the gamma event time intervals (iv) for
%     one band (ie evts.low)
%
%   Outputs:
%     - data_ft [struct]: data split into trials in the field trip format.
%
%   Requires FieldTrip toolbox for "ft_redefinetrial"
%
%
%% Defaults

%cfg_def.buffer = 0.25; % buffer of 250ms on either side of the center of the event. Used for plots
cfg_def.twin = 0.025;  % default for low gamma
cfg_def.check = 0;
cfg = ProcessConfig(cfg_def, cfg_in);

%% trialify
event_centers = IVcenters(evts_in);

cfg_trl = [];
cfg_trl.t = cat(1,event_centers);
cfg_trl.t = cfg_trl.t - evts_in.firstTimestamp;
cfg_trl.twin = [-cfg.twin cfg.twin];
cfg_trl.hdr = data_remap_FT.hdr;

trl = ft_maketrl(cfg_trl);

cfg_trl = [];
cfg_trl.trl = trl;
data_trl = ft_redefinetrial(cfg_trl,data_remap_FT);

%% Check the data

if cfg.check ==1
    AMPX_trial_check(data_trl)
end