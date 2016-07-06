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
extract_varargin

%%
% CD to the data directory.  Should contain all the .dat, .ini, and .meta
fname = strrep(fname, '_', '-');
cd(['D:\DATA\' fname(1:4) '\' fname(1:15) ])
% mkpre % gets the file name for the pre record.
%%  determine session type
% is this sessions a pre, task or post record?
% session_type = type;

if isempty(str2num(fname(end)))==1
    fname = fname(1:15);
end

%% preprocess/load the data

[data, data_ft] = AMPX_Naris_preprocess([],fname,session_type);

LoadExpKeys
%% Organize the data to fit the probe lay out in the 8x8 format.
% Data about the probemapping

data_remap_AMPX = AMPX_remap_sites(data, ExpKeys)

data_remap_FT = AMPX_remap_sites(data_ft, ExpKeys)
clear data_ft

ExpKeys_remap = AMPX_remap_ExpKeys(ExpKeys, data_remap_AMPX);
%% get the events (using Catanese et al. method)
if strcmp(fname(1:4), 'R045') % this channel is cleaner for detection but not the most lateral thus is it different than what is used for phase analysis later on. 
[evts, ~, detect_chan] = AMPX_Julien_DetectEvents(data_remap_AMPX, ExpKeys_remap, 'debug2', 1, 'detect_chan', 48);
else
    [evts, ~, detect_chan] = AMPX_Julien_DetectEvents(data_remap_AMPX, ExpKeys_remap, 'debug2', 1);
end

% for R045 use the 98th percentile instead of the 95th for high gamma

if strcmp(fname(1:4), 'R045')
%     evts.low = evts.low_98;
    evts.high = evts.high_98;
end

if strcmp(session_type, 'task')
    feeder_evts = AMPX_get_task_events([], data_remap_AMPX);
    evts.feeder.low = IntersectIV([], evts.low,feeder_evts.feeder_iv);
    evts.feeder.high = IntersectIV([], evts.high, feeder_evts.feeder_iv);
    evts.approach.low = IntersectIV([],  evts.low, feeder_evts.approach_iv);
    evts.approach.high = IntersectIV([],  evts.high, feeder_evts.approach_iv);
    evts.photo.low = IntersectIV([],  evts.low, feeder_evts.photo_iv);
    evts.photo.high = IntersectIV([], evts.high,feeder_evts.photo_iv);
    
    fprintf(['Feeder Low: ' num2str(length(evts.feeder.low.tstart)) '   High: ' num2str(length(evts.feeder.high.tstart)) '\n'])
    fprintf(['Approach Low: ' num2str(length(evts.approach.low.tstart)) '   High: ' num2str(length(evts.approach.high.tstart)) '\n'])
    fprintf(['Photo Low: ' num2str(length(evts.photo.low.tstart)) '   High: ' num2str(length(evts.photo.high.tstart)) '\n'])
end

%% get an equal number of non-overlapping
match_tsd = AMPX_to_tsd(data_remap_AMPX);
match_tsd.data = match_tsd.data(detect_chan, :);
match_tsd.label = detect_chan;
cfg = []; cfg.debug = 1; cfg.PARAM_control_dt = 20;
ctrl_evts = AMPX_match_events(cfg, evts, match_tsd);
close all
%% consolidate all the events
evts.ctrl_low = ctrl_evts.lg;
evts.ctrl_high = ctrl_evts.hg;

evts.ctrl_low.firstTimestamp = evts.low.firstTimestamp;
evts.ctrl_high.firstTimestamp = evts.high.firstTimestamp;

%% check
% match_tsd = AMPX_to_tsd(data);
% match_tsd.data = match_tsd.data(67, :);
% match_tsd.label = 67;
% cfg_plot = []; cfg_plot.display = 'iv'; cfg_plot.mode = 'center';
% PlotTSDfromIV(cfg_plot, evts.high, match_tsd)

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

if strcmp(session_type, 'task')
    % low reward
    cfg = [];
    cfg.twin = 0.05;
    cfg.check = 1; % make a plot of 16 random trials to see the quality
    data_low_reward_ft = AMPX_to_FT_trial_split(cfg, evts.feeder.low, data_remap_FT);
    
    % low approach
    cfg = [];
    cfg.twin = 0.05;
    cfg.check = 1; % make a plot of 16 random trials to see the quality
    data_low_approach_ft = AMPX_to_FT_trial_split(cfg, evts.approach.low, data_remap_FT);
    
    % high reward
    cfg = [];
    cfg.twin = 0.05;
    cfg.check = 1; % make a plot of 16 random trials to see the quality
    data_high_reward_ft = AMPX_to_FT_trial_split(cfg, evts.feeder.high, data_remap_FT);
    
    % high approach
    cfg = [];
    cfg.twin = 0.05;
    cfg.check = 1; % make a plot of 16 random trials to see the quality
    data_high_approach_ft = AMPX_to_FT_trial_split(cfg, evts.approach.high, data_remap_FT);
end

% Spindles
if isempty(evts.spindles.tstart) ~=1
    cfg = [];
    cfg.twin = 0.25;
    cfg.check = 1; % make a plot of 9 random trials to see the quality
    data_spindles_ft = AMPX_to_FT_trial_split(cfg, evts.spindles, data_remap_FT);
end
close all

%% get power

% low gamma
cfg = [];
cfg.freq = [40 55];
data_out.lg.power = AMPX_get_power(cfg, data_lg_ft, ExpKeys_remap);

% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg.power = AMPX_get_power(cfg, data_hg_ft, ExpKeys_remap);

% low random
cfg = [];
cfg.freq = [40 55];
data_out.lg_ran.power = AMPX_get_power(cfg, data_ctrl_lg_ft, ExpKeys_remap);

% high gamma
cfg = [];
cfg.freq = [70 85];
data_out.hg_ran.power = AMPX_get_power(cfg, data_ctrl_hg_ft, ExpKeys_remap);

if strcmp(session_type, 'task')
    % low gamma
    if isempty(evts.feeder.low.tstart) ==1
        data_out.lg_reward.power.power = [];
        data_out.lg_reward.power.power_avg = [];
    else
        cfg = [];
        cfg.freq = [40 55];
        data_out.lg_reward.power = AMPX_get_power(cfg, data_low_reward_ft, ExpKeys_remap);
    end
    % high gamma
    if isempty(evts.feeder.high.tstart) ==1
        data_out.hg_reward.power.power = [];
        data_out.hg_reward.power.power_avg = [];
    else
        cfg = [];
        cfg.freq = [70 85];
        data_out.hg_reward.power = AMPX_get_power(cfg, data_high_reward_ft, ExpKeys_remap);
    end
    % low approach
    if isempty(evts.approach.low.tstart) ==1
        data_out.lg_approach.power.power = [];
        data_out.lg_approach.power.power_avg = [];
    else
        cfg = [];
        cfg.freq = [40 55];
        data_out.lg_approach.power = AMPX_get_power(cfg, data_low_approach_ft, ExpKeys_remap);
    end
    % high gamma
    if isempty(evts.approach.low.tstart) ==1
        data_out.hg_approach.power.power = [];
        data_out.hg_approach.power.power_avg = [];
    else
        cfg = [];
        cfg.freq = [70 85];
        data_out.hg_approach.power = AMPX_get_power(cfg, data_high_approach_ft, ExpKeys_remap);
    end
end
%% get the plane fitting statistics
if strcmp(session_type, 'task') ~=1
    cfg= [];
    cfg.name = fname;
    cfg.session_type = session_type;
    
    % low gamma
    disp('low gamma')
    cfg.type = 'low';
    data_out.lg.power.stats = AMPX_get_plane_fitting_advanced(cfg,data_out.lg.power);
    plane_stats.low_gamma.rsq =data_out.lg.power.stats.rsq;
    plane_stats.low_gamma.rsq2 =data_out.lg.power.stats.rsq2;
    
    % high gamma
    disp('high gamma')
    cfg.type = 'high';
    data_out.hg.power.stats = AMPX_get_plane_fitting_advanced(cfg, data_out.hg.power);
    plane_stats.high_gamma.rsq =data_out.hg.power.stats.rsq;
    plane_stats.high_gamma.rsq2 =data_out.hg.power.stats.rsq2;
    
    % low control
    disp('low control')
    cfg.type = 'low_ctrl';
    data_out.lg_ran.power.stats = AMPX_get_plane_fitting_advanced(cfg, data_out.lg_ran.power);
    plane_stats.low_control.rsq =data_out.lg_ran.power.stats.rsq;
    plane_stats.low_control.rsq2 =data_out.lg_ran.power.stats.rsq2;
    
    % high control
    disp('high control')
    cfg.type = 'high_ctrl';
    data_out.hg_ran.power.stats = AMPX_get_plane_fitting_advanced(cfg, data_out.hg_ran.power);
    plane_stats.high_control.rsq =data_out.hg_ran.power.stats.rsq;
    plane_stats.high_control.rsq2 =data_out.hg_ran.power.stats.rsq2;
    
    
    %% smoothed plane fitting
    cfg = [];
    cfg.name = fname;
    cfg.session_type = session_type;
    cfg.plot = 1;
    % low gamma
    disp('low gamma')
    cfg.type = 'low';
    data_out.lg.power.smooth.stats = AMPX_get_plane_fitting(cfg,data_out.lg.power);
    plane_stats.low_gamma.smooth.rsq =data_out.lg.power.smooth.stats.rsq;
    plane_stats.low_gamma.smooth.rsq2 =data_out.lg.power.smooth.stats.rsq2;
    
    % high gamma
    disp('high gamma')
    cfg.type = 'high';
    data_out.hg.power.smooth.stats = AMPX_get_plane_fitting(cfg, data_out.hg.power);
    plane_stats.high_gamma.smooth.rsq =data_out.hg.power.smooth.stats.rsq;
    plane_stats.high_gamma.smooth.rsq2 =data_out.hg.power.smooth.stats.rsq2;
    
    % low control
    disp('low control')
    cfg.type = 'low_ctrl';
    data_out.lg_ran.power.smooth.stats = AMPX_get_plane_fitting(cfg, data_out.lg_ran.power);
    plane_stats.low_control.smooth.rsq =data_out.lg_ran.power.smooth.stats.rsq;
    plane_stats.low_control.smooth.rsq2 =data_out.lg_ran.power.smooth.stats.rsq2;
    
    % high control
    disp('high control')
    cfg.type = 'high_ctrl';
    data_out.hg_ran.power.smooth.stats = AMPX_get_plane_fitting(cfg, data_out.hg_ran.power);
    plane_stats.high_control.smooth.rsq =data_out.hg_ran.power.smooth.stats.rsq;
    plane_stats.high_control.smooth.rsq2 =data_out.hg_ran.power.smooth.stats.rsq2;
    
    %% get the phase differences across an entire event
    % low gamma
    cfg =[];
    cfg.freq = [40 55];
    data_out.lg.phase = AMPX_get_phase(cfg, data_lg_ft, ExpKeys_remap);
    
    
    % high gamma
    cfg =[];
    cfg.freq = [70 85];
    data_out.hg.phase = AMPX_get_phase(cfg, data_hg_ft, ExpKeys_remap);
    
    % low control
    cfg =[];
    cfg.freq = [40 55];
    data_out.lg_ran.phase = AMPX_get_phase(cfg, data_ctrl_lg_ft, ExpKeys_remap);
    
    
    % high control
    cfg =[];
    cfg.freq = [70 85];
    data_out.hg_ran.phase = AMPX_get_phase(cfg, data_ctrl_hg_ft, ExpKeys_remap);
    
    %% extract middle three cycles in each gamma event
    
    % low gamma
    cfg = []; cfg.freq = [40 55];
    data_out.lg.cycles = AMPX_get_3cycles(cfg, evts.low, data_remap_AMPX, ExpKeys_remap);
    
    
    % high gamma
    cfg =[];
    cfg.freq = [70 85];
    data_out.hg.cycles = AMPX_get_3cycles(cfg, evts.high, data_remap_AMPX, ExpKeys_remap);
    
    % low control
    cfg = [];
    cfg.freq = [40 55];
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
    
    % %% Prepare the data_out struct for kCSD
    % cfg = [];
    % data_out.lg.kCSD = AMPX_cycle_kCSD(cfg, data_out.lg.cycles);
    %% put the events back in the output structure
    % low
    data_out.lg.evts = evts.low;
    
    %high
    data_out.hg.evts = evts.high;
    
    % low control
    data_out.lg_ran.evts = evts.ctrl_low;
    
    % high control
    data_out.hg_ran.evts = evts.ctrl_high;
    
    % spindles
    data_out.spindles.evts = evts.spindles;
    
end
if strcmp(session_type, 'task')
    % low feeder
    data_out.lg_reward.evts = evts.feeder.low;
    % low approach
    data_out.lg_approach.evts = evts.approach.low;
    %high feeder
    data_out.hg_reward.evts = evts.feeder.high;
    % high approach
    data_out.hg_approach.evts = evts.approach.high;
end

% save(['Naris_data_out_' session_type], 'data_out')
