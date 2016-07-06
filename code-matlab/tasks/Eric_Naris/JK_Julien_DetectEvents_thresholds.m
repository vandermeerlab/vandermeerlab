function [evts_out, csc, evt_threshold] = JK_Julien_DetectEvents_thresholds(cfg_in, data, ExpKeys)

%% MASTER_CollectGammaEvents.m
% detects and organizes gamma events from raw LFP data
%
% Julien Catanese & Matthijs van der Meer

%% Convert data to TSD and extract the channel for detection from the ExpKeys

% data_tsd = AMPX_to_tsd(data);
csc = data;
% csc.data = csc.data(ExpKeys.DetectChan,:);
% csc.detect_chan = data% ExpKeys.DetectChan;
%% set params
% gamma event detection
cfg_def.f_label = {'low','high', 'low_low_tr', 'high_low_tr'};
cfg_def.f_bandpass = {[40 55],[70 85],[40 55], [70 85]}; % frequency bands for event detection
cfg_def.detect_thr = [0.95 .95 0.90 .90]; % threshold for event detection: 95th percentile of (amplitude) envelope
cfg_def.detect_method = 'percentile'; %'raw'; % 'raw', 'zscore', 'percentile'
cfg_def.detect_nCycles = 4; % require minimum number of gamma cycles
cfg_def.var_thr = 1.5; % def: variance/mean of cycle peaks and throughs must be smaller than this
cfg_def.detect_epoch = 'all'; % 'all', 'post', 'task'; % set threshold based on what data (for events)
cfg_def.ampl_min = [0.25e-04 0.2e-04 0.25e-04 0.2e-04]; % all peaks of event must be larger than this (in V)

% artifact, chewing, and spindle detection
cfg_def.artif_thr =  std(csc.data)*4;   %0.75 * 10^-3; % raw amplitude must be smaller than this (in V) to pass artifact detection
cfg_def.chew_thr = 3; % z-score of chew band envelope must be smaller than this, default 0.25
cfg_def.spindle_thr = 4; % z-score of spindle band envelope must be smaller than this

cfg = ProcessConfig2(cfg_def, cfg_in);
% some flags to enable visualization
debug = 0; debug2 = 0;

%% main loop over sessions


%% detect major transient artifacts: this is on UNFILTERED data because of sustained artifacts in R026 (see project log)
csc_artif = csc;
csc_artif.data = abs(csc_artif.data); % detect artifacts both ways

cfg_artif_det = [];
cfg_artif_det.method = 'raw';
cfg_artif_det.threshold = cfg.artif_thr;
cfg_artif_det.minlen = 0;
evt_artif = TSDtoIV(cfg_artif_det,csc_artif);

cfg_temp = []; cfg_temp.d = [-0.5 0.5];
evt_artif = ResizeIV(cfg_temp,evt_artif);

%% plot
if debug
    cfg_plot.display = 'tsd'; % tsd, iv
    PlotTSDfromIV(cfg_plot,evt_artif,csc);
    pause; close all;
end

%% NaN out artifacts to improve reliability of subsequent z-scoring
artif_idx = TSD_getidx2(csc,evt_artif); % if error, try TSD_getidx (slower)
csc.data(artif_idx) = NaN;

%% find chewing artifacts
cfg_chew = [];
cfg_chew.epoch = 'all'; % chewing occurs during task mostly
cfg_chew.minlen = 0.02;
cfg_chew.filter_cfg.f = [200 300]; % default [200 300]
cfg_chew.threshold = cfg.chew_thr; % 0.25 for session 2, 0.5 for session 1?
cfg_chew.smooth = 0.05; % convolve with Gaussian of this SD
evt_chew = Julien_DetectEvents(cfg_chew,csc,ExpKeys);

cfg_temp = []; cfg_temp.d = [-0.1 0.1];
evt_chew = ResizeIV(cfg_temp,evt_chew);

%% plot
if debug
    cfg_plot.display = 'tsd'; % tsd, iv
    PlotTSDfromIV(cfg_plot,evt_chew,csc);
    pause; close all;
end

%% NaN out artifacts to improve reliability of subsequent z-scoring
chew_idx = TSD_getidx2(csc,evt_chew); % if error, try TSD_getidx (slower)
csc.data(chew_idx) = NaN;
%% find sleep spindles
cfg_spindl = [];
cfg_spindl.epoch = 'all'; % chewing occurs during task mostly
cfg_spindl.minlen = 0.005;
cfg_spindl.filter_cfg.f = [25 30]; % default [25 35]
cfg_spindl.threshold = cfg.spindle_thr; % 0.25 for session 2, 0.5 for session 1?
cfg_spindl.smooth = 0.05; % convolve with Gaussian of this SD
evt_spindl = Julien_DetectEvents(cfg_spindl,csc,ExpKeys);

cfg_temp = []; cfg_temp.d = [-0.5 0.5];
evt_spindl = ResizeIV(cfg_temp,evt_spindl);

%% plot
if debug
    cfg_plot.display = 'tsd'; % tsd, iv
    PlotTSDfromIV(cfg_plot,evt_spindl,csc);
    pause; close all;
end

%% now, loop over frequency bands to process
for iFreq = 1:length(cfg.f_label)
    
    %% set up filter
    cfg_filter = [];
    cfg_filter.f = cfg.f_bandpass{iFreq};
    cfg_filter.type = 'cheby1';
    cfg_filter.order = 5;
    
    %% basic gamma detection
    cfg_evt = [];
    cfg_evt.epoch = cfg.detect_epoch;
    cfg_evt.epochLength = 10*60; % was 5 * 60
    cfg_evt.filter_cfg = cfg_filter;
    cfg_evt.minlen = cfg.detect_nCycles./mean(cfg.f_bandpass{iFreq}); % or, 0.05
    %cfg_evt.minlen = 0.05;
    cfg_evt.smooth = 0.05; % convolve with Gaussian of this SD
    cfg_evt.threshold = cfg.detect_thr(iFreq);
    cfg_evt.method = cfg.detect_method;
    
    [evts,evt_thr] = Julien_DetectEvents(cfg_evt,csc,ExpKeys);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events detected initially.\n',length(evts.tstart),cfg.f_label{iFreq});
    
    if debug
        cfg_plot = []; cfg_plot.display = 'iv'; cfg_plot.mode = 'center'; cfg_plot.width = 0.2;
        PlotTSDfromIV(cfg_plot,evts,csc);
        pause; close all;
    end
   
   evt_threshold(iFreq) = evt_thr; 
   evts_out.(cfg.f_label{iFreq}) = evts;
%    detect_csc = [] ;
end