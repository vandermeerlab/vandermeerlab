function [data_out] = AMPX_Naris_pipeline(fname, varargin)
%% AMPX_Naris_pipeline performs: powerspectral densities, phase offset and PCA on a given data session.
%                          This uses:
%                          AMPX_loadData
%                          AMPX_RemoveArtifacts
%                          AMPX_filter
%                          AMPX_trial_split
%                          AMPX_get_power
%                          AMPX_get_phase
%                          AMPX_phase_distrib_plot
%
%
%INPUTS:
%   - file_name: 'R0XX-yyyy-mm-dd' If this is a pre or post session add
%   "-post" or "-pre"
% Optional:
%   - 'gamma_type' : either 'high', 'low' or if left empty will do both
%   - 'figures' : 'on' or 'off' [default: 'off']
%   - 'savefig' : 'on' or 'off' [default: 'off']
%   - 'gamma_chan' : channel for gamma detection.  This is normally done
%   using the channel with the highest power.
%OUTPUTS:
%   - 'data_out' : a structure containing the power/phase/PCA results for
%   both the high and low gamma events as well as some of the parameters
%   used.
%
% EC 04/2016


%% Preamble

% CD to the data directory.  Should contain all the .dat, .ini, and .meta
cd('D:\DATA\R054\R054-2014-10-10')
mkpre % gets the file name for the pre record.  
%%  determine session type
fname = strrep(file_name, '_', '-');
% is this sessions a pre, task or post record?
if strcmp(fname(end-3:end), 'post')
    session_type = 'post';
elseif strcmp(fname(end-2:end), 'pre')
    session_type = 'pre';
else
    session_type = 'task';
end

if isempty(str2num(fname(end)))==1
    fname = fname(1:15);
end

%% preprocess/load the data

[data, data_ft] = AMPX_Naris_preprocess([],file_name,session_type);

LoadExpKeys

%% get the events (using Catanese et al. method)

[evts, ~, detect_chan] = AMPX_Julien_DetectEvents(data, ExpKeys); 

%% get an equal number of non-overlapping 
match_tsd = AMPX_to_tsd(data);
match_tsd.data = match_tsd.data(detect_chan, :);
match_tsd.label = detect_chan;
cfg = []; cfg.debug = 1; cfg.PARAM_control_dt = 10;
ctrl_evts = AMPX_match_events(cfg, evts, match_tsd);

%% consolidate all the events
evts.ctrl_low = ctrl_evts.lg;
evts.ctrl_high = ctrl_evts.hg;

evts.ctrl_low.firstTimestamp = evts.low.firstTimestamp;
evts.ctrl_high.firstTimestamp = evts.high.firstTimestamp;

%% check 
cfg_plot = []; cfg_plot.display = 'iv'; cfg_plot.mode = 'center';
PlotTSDfromIV(cfg_plot, evts.ctrl_low, match_tsd)
%% Organize the data to fit the probe lay out in the 8x8 format.
% Data about the probemapping 

data_remap_AMPX = AMPX_remap_sites(data, ExpKeys)

data_remap_FT = AMPX_remap_sites(data_ft, ExpKeys)

ExpKeys_remap = AMPX_remap_ExpKeys(ExpKeys, data_remap_AMPX);

%% define FT trials based on evts
addpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip')
ft_defaults
% low gamma
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_lg_ft = AMPX_to_FT_trial_split(cfg, evts.low, data_remap_FT);


% high Gamma
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_hg_ft = AMPX_to_FT_trial_split(cfg, evts.high, data_remap_FT);


% low control
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_ctrl_lg_ft = AMPX_to_FT_trial_split(cfg, evts.ctrl_low, data_remap_FT);


% high Control
cfg = [];
cfg.twin = 0.05;
cfg.check = 1; % make a plot of 16 random trials to see the quality
data_ctrl_hg_ft = AMPX_to_FT_trial_split(cfg, evts.ctrl_high, data_remap_FT);


% Spindles
cfg = [];
cfg.twin = 0.25;
cfg.check = 1; % make a plot of 9 random trials to see the quality
data_spindles_ft = AMPX_to_FT_trial_split(cfg, evts.spindles, data_remap_FT);

%% get power

% low gamma
cfg = [];
cfg.freq = [40 55];
data_out.lg.power = AMPX_get_power(cfg, evts.low, data_lg_ft, ExpKeys_remap);

% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg.power = AMPX_get_power(cfg, evts.high, data_hg_ft, ExpKeys_remap);

% low random
cfg = [];
cfg.freq = [40 55];
data_out.lg_ran.power = AMPX_get_power(cfg, evts.ctrl_low, data_ctrl_lg_ft, ExpKeys_remap);

% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg_ran.power = AMPX_get_power(cfg, evts.ctrl_high, data_ctrl_hg_ft, ExpKeys_remap);

%% get the plane fitting statistics 
% low gamma
disp('low gamma')
data_out.lg.power.stats = AMPX_get_plane_fitting(data_out.lg.power);
plane_stats.low_gamma.rsq =data_out.lg.power.stats.rsq;
plane_stats.low_gamma.rsq2 =data_out.lg.power.stats.rsq2;

% high gamma
disp('high gamma')
data_out.hg.power.stats = AMPX_get_plane_fitting(data_out.hg.power);
plane_stats.high_gamma.rsq =data_out.hg.power.stats.rsq;
plane_stats.high_gamma.rsq2 =data_out.hg.power.stats.rsq2;

% low control
disp('low control')
data_out.lg_ran.power.stats = AMPX_get_plane_fitting(data_out.lg_ran.power);
plane_stats.low_control.rsq =data_out.lg_ran.power.stats.rsq;
plane_stats.low_control.rsq2 =data_out.lg_ran.power.stats.rsq2;

% high control
disp('high control')
data_out.hg_ran.power.stats = AMPX_get_plane_fitting(data_out.hg_ran.power);
plane_stats.high_control.rsq =data_out.hg_ran.power.stats.rsq;
plane_stats.high_control.rsq2 =data_out.hg_ran.power.stats.rsq2;

% output the plane stats
AMPX_plane_count_hist(plane_stats.low_gamma, plane_stats.low_control, 'low')

AMPX_plane_count_hist(plane_stats.high_gamma, plane_stats.high_control, 'high')

%% get the phase differences across an entire event
% low gamma
cfg =[];
cfg.freq = [40 55];
data_out.lg.phase = AMPX_get_phase(cfg, data_lg_ft, ExpKeys_remap);


% high gamma
cfg =[];
cfg.freq = [70 85];
data_out.hg.power = AMPX_get_phase(cfg, data_hg_ft, ExpKeys_remap); 

% low control
cfg =[];
cfg.freq = [40 55];
data_out.lg_ran.phase = AMPX_get_phase(cfg, data_ctrl_lg_ft, ExpKeys_remap);


% high control
cfg =[];
cfg.freq = [70 85];
data_out.hg_ran.power = AMPX_get_phase(cfg, data_ctrl_hg_ft, ExpKeys_remap);

%% extract middle three cycles in each gamma event

% low gamma
cfg = []; cfg.freq = [40 55];
data_out.lg.cycles = AMPX_get_3cycles(cfg, evts.low, data_remap_AMPX, ExpKeys_remap);


% high gamma
cfg =[];
cfg.freq = [70 85];
data_out.hg.cycles = AMPX_get_3cycles(cfg, evts.high, data_remap_AMPX, ExpKeys_remap);

% low control
cfg = []; cfg.freq = [40 55];

data_out.lg_ran.cycles = AMPX_get_3cycles(cfg, evts.low, data_remap_AMPX, ExpKeys_remap);


% high control
cfg =[];
cfg.freq = [70 85];
data_out.hg_ran.cycles = AMPX_get_3cycles(cfg, evts.ctrl_high, data_remap_AMPX, ExpKeys_remap);

%% get the phase for the cycle data
% low gamma
cfg = []; cfg.f= [40 55];
data_out.lg.cycles_phase = AMPX_phase_cycle(cfg, data_out.lg.cycles);

% high gamma
cfg =[];
cfg.f= [70 85];
data_out.hg.cycles_phase = AMPX_phase_cycle(cfg, data_out.hg.cycles);
% low control 
cfg = []; cfg.f = [40 55];
data_out.lg_ran.cycles_phase = AMPX_phase_cycle(cfg, data_out.lg_ran.cycles);

% high control
cfg =[];
cfg.f = [70 85];
data_out.hg_ran.cycles_phase = AMPX_phase_cycle(cfg, data_out.hg_ran.cycles);

%% Prepare the data_out struct for kCSD
cfg = [];
data_out.lg.kCSD = AMPX_cycle_kCSD(cfg, data_out.lg.cycles.peaks);

