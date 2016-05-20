function [evts, detect_csc, detect_chan] = AMPX_Julien_DetectEvents(data, ExpKeys)

%% MASTER_CollectGammaEvents.m
% detects and organizes gamma events from raw LFP data
%
% Julien Catanese & Matthijs van der Meer

%% Convert data to TSD and extract the channel for detection from the ExpKeys

data_tsd = AMPX_to_tsd(data);
csc = data_tsd;
csc.data = csc.data(ExpKeys.DetectChan,:);
csc.detect_chan = ExpKeys.DetectChan;
detect_chan = ExpKeys.DetectChan;
%% set params
% gamma event detection
PARAM_f_label = {'low','high', 'low_low_tr', 'high_low_tr'};
PARAM_f_bandpass = {[40 55],[70 85],[40 55], [70 85]}; % frequency bands for event detection
PARAM_detect_thr = [0.95 .95 0.90 .90]; % threshold for event detection: 95th percentile of (amplitude) envelope
PARAM_detect_method = 'percentile'; % 'raw', 'zscore', 'percentile'
PARAM_detect_nCycles = 4; % require minimum number of gamma cycles
PARAM_var_thr = 1.5; % def: variance/mean of cycle peaks and throughs must be smaller than this
PARAM_detect_epoch = 'all'; % 'all', 'post', 'task'; % set threshold based on what data (for events)
PARAM_ampl_min = [0.5e-04 0.5e-04 0.5e-04 0.5e-04]; % all peaks of event must be larger than this (in V)

% artifact, chewing, and spindle detection
PARAM_artif_thr =  std(csc.data)*4;   %0.75 * 10^-3; % raw amplitude must be smaller than this (in V) to pass artifact detection
PARAM_chew_thr = 3; % z-score of chew band envelope must be smaller than this, default 0.25
PARAM_spindle_thr = 4; % z-score of spindle band envelope must be smaller than this

% some flags to enable visualization
debug = 0; debug2 = 0;

%% main loop over sessions


%% detect major transient artifacts: this is on UNFILTERED data because of sustained artifacts in R026 (see project log)
csc_artif = csc;
csc_artif.data = abs(csc_artif.data); % detect artifacts both ways

cfg_artif_det = [];
cfg_artif_det.method = 'raw';
cfg_artif_det.threshold = PARAM_artif_thr;
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
cfg_chew.filter_cfg.f = [200 500]; % default [200 300]
cfg_chew.threshold = PARAM_chew_thr; % 0.25 for session 2, 0.5 for session 1?
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
cfg_spindl.threshold = PARAM_spindle_thr; % 0.25 for session 2, 0.5 for session 1?
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
for iFreq = 1:length(PARAM_f_label)
    
    %% set up filter
    cfg_filter = [];
    cfg_filter.f = PARAM_f_bandpass{iFreq};
    cfg_filter.type = 'cheby1';
    cfg_filter.order = 5;
    
    %% basic gamma detection
    cfg_evt = [];
    cfg_evt.epoch = PARAM_detect_epoch;
    cfg_evt.epochLength = 10*60; % was 5 * 60
    cfg_evt.filter_cfg = cfg_filter;
    cfg_evt.minlen = PARAM_detect_nCycles./mean(PARAM_f_bandpass{iFreq}); % or, 0.05
    %cfg_evt.minlen = 0.05;
    cfg_evt.smooth = 0.05; % convolve with Gaussian of this SD
    cfg_evt.threshold = PARAM_detect_thr(iFreq);
    cfg_evt.method = PARAM_detect_method;
    
    [evt,evt_thr] = Julien_DetectEvents(cfg_evt,csc,ExpKeys);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events detected initially.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    if debug
        cfg_plot = []; cfg_plot.display = 'iv'; cfg_plot.mode = 'center'; cfg_plot.width = 0.2;
        PlotTSDfromIV(cfg_plot,evt,csc);
        pause; close all;
    end
    
    %% remove artifacts
    evt = DifferenceIV([],evt,evt_artif);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after artifact removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    %% remove chewing
    evt = DifferenceIV([],evt,evt_chew);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after chewing removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    %% remove spindles
    evt = DifferenceIV([],evt,evt_spindl);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after spindle removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    
    %% exclude events with insufficient gamma cycles - count how many exist above same threshold as used for detection
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_filter;
    evt = CountCycles(cfg_cc,csc,evt);
    
    cfg_cc = [];
    cfg_cc.operation = '>=';
    cfg_cc.threshold = PARAM_detect_nCycles-1;
    evt = SelectIV(cfg_cc,evt,'nCycles');
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after cycle count removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    %% exclude events with excessive variance in amplitude (only those in within evt detect boundaries)
    %         iv_temp = IVcenters(evt);
    %         iv_temp = iv(iv_temp-cfg_evt.minlen/2,iv_temp+cfg_evt.minlen/2);
    %
    %         cfg_cc = [];
    %         cfg_cc.threshold_type = 'raw';
    %         cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    %         cfg_cc.filter_cfg = cfg_filter;
    %         iv_temp = CountCycles(cfg_cc,csc,iv_temp);
    %
    %         cfg_cc = [];
    %         cfg_cc.operation = '<';
    %         cfg_cc.threshold = PARAM_var_thr;
    %         [~,keep_idx] = SelectIV(cfg_cc,iv_temp,'var_raw');
    %
    %         evt = SelectIV([],evt,keep_idx); iv_temp = SelectIV([],iv_temp,keep_idx);
    %
    %         cfg_cc = [];
    %         cfg_cc.operation = '<';
    %         cfg_cc.threshold = PARAM_var_thr;
    %         [~,keep_idx] = SelectIV(cfg_cc,iv_temp,'var');
    %
    %         evt = SelectIV([],evt,keep_idx);
    
    
    %% exclude events with excessive variance in amplitude (all peaks and troughs)
    cfg_cc = [];
    cfg_cc.operation = '<';
    cfg_cc.threshold = PARAM_var_thr;
    evt = SelectIV(cfg_cc,evt,'var_raw');
    
    cfg_cc = [];
    cfg_cc.operation = '<';
    cfg_cc.threshold = PARAM_var_thr;
    evt = SelectIV(cfg_cc,evt,'var');
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after amplitude variance removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    %% exclude events with insufficient mean (or min) -- could try doing this on minlen part only
    %         cfg_cc = [];
    %         cfg_cc.operation = '>';
    %         cfg_cc.threshold = PARAM_ampl_min(iFreq);
    %         evt = SelectIV(cfg_cc,evt,'min_filt');
    %
    %         fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after mean peak removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    %
    %% minlen only version
    iv_temp = IVcenters(evt);
    iv_temp = iv(iv_temp-cfg_evt.minlen/2,iv_temp+cfg_evt.minlen/2);
    
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_filter;
    iv_temp = CountCycles(cfg_cc,csc,iv_temp);
    
    cfg_cc = [];
    cfg_cc.operation = '>';
    cfg_cc.threshold = PARAM_ampl_min(iFreq);
    [~,keep_idx] = SelectIV(cfg_cc,iv_temp,'min_filt');
    
    evt = SelectIV([],evt,keep_idx);
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after mean peak removal.\n',length(evt.tstart),PARAM_f_label{iFreq});
    
    
    %% visualize
    if debug2
        % raw LFP only
        cfg_plot = []; cfg_plot.display = 'iv'; cfg_plot.mode = 'center'; cfg_plot.width = 0.2; cfg_plot.title = 'var_raw';
        PlotTSDfromIV(cfg_plot,evt,csc);
        pause; close all;
        
        % TFR version
        %             cfg_convert = []; cfg_convert.mode = 'resample';
        %             csc_ft = TSDtoFT(cfg_convert,csc);
        %
        %             evt_temp = evt;
        %             evt_temp.tstart = evt_temp.tstart-csc.tvec(1);
        %             evt_temp.tend = evt_temp.tend-csc.tvec(1);
        %
        %             cfg_temp = []; cfg_temp.foi = 1:5:300; cfg_temp.clim = [0 10^-9]; cfg_temp.twin = [-0.1 0.1];
        %             PlotTSDfromIV_TFR(cfg_temp,evt_temp,csc_ft);
        %             pause; close all;
    end
    
    %% store data into collector variable
    evts.(PARAM_f_label{iFreq}) = evt;
    evts.(PARAM_f_label{iFreq}).firstTimestamp = csc.tvec(1);
    
    %         sess_id_field = regexprep(SessionLIST{iFD},'-','_');
    %         ALL_evt.(this_RatID).(sess_id_field).(PARAM_f_label{iFreq}) = evt;
    %         ALL_evt.(this_RatID).(sess_id_field).fd = fd{iFD};
    %         ALL_evt.(this_RatID).(sess_id_field).firstTimestamp = csc.tvec(1); % need this for trialification later
    
end % of loop over frequencies

%% remove overlapping gamma events.  
temp_lg = []; temp_hg = []; temp_lg_low = []; temp_hg_low = [];

temp_lg = DifferenceIV([],evts.low,evts.high);
temp_hg = DifferenceIV([],evts.high,evts.low);

temp_lg_low = DifferenceIV([],evts.low_low_tr,evts.high_low_tr);
temp_hg_low = DifferenceIV([],evts.high_low_tr,evts.low_low_tr);


evts.low = temp_lg;
evts.high = temp_hg;
evts.high_tr = temp_hg_low;
evts.low_tr = temp_lg_low;
evts.spindles = evt_spindl;
evts.spindles.firstTimestamp = csc.tvec(1);
%%


for iFreq = 1:length(PARAM_f_label)
    fprintf(['\n ' PARAM_f_label{iFreq} ' found: ' num2str(length(evts.(PARAM_f_label{iFreq}).tstart)) ' events'])
end

fprintf(['\nSpindles found:'  num2str(length(evts.spindles.tstart)) ' events'])
fprintf('\n')

% output the same csc used to detect the events (used for pseudo random
% non-gamma epochs selection
detect_csc = csc;
