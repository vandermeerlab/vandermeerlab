function spwr_iv_out = getSWR(cfg_in,S,lfp)
%% GETSWR Get Sharp Wave Ripple events
%   spwr_iv_out = getSWR(cfg_in,S,lfp) returns an IV struct containing the time intervals
%   of sharp wave ripple events found in lfp. GETSWR uses both lfp power thresholding and
%   multi-unit activity to determine valid time intervals.
%
%   INPUT:
%       cfg - input cfg parameters
%       S - input spiketrain
%       lfp - input csc
%
%   OUTPUT:
%       swr_iv - iv containing swr events
%
%   CFG OPTIONS:
%       cfg.f = [140 200]; % passband range in Hz
%       cfg.threshold = 3; %in std, determines size of ivs
%       cfg.max_thr = 5; % in std, for selecting ivs
%       cfg.merge_thr = 0.03;
%       cfg.viewMUA = 0; % display MUA
%       cfg.useMUA = 0; % use MUA for SWR detection
%       cfg.overlap_thr = 0.01; %in seconds, overlap between MUA and SWR
%
% youkitan 2014-12-28

%% Parse cfg parameters
cfg_def.f = [140 200]; % ripple-band power passband (Hz)
cfg_def.threshold = 2; % iv edge threshold for ripple-band power (z-score)
cfg_def.max_thr = 5; % iv max threshold for ripple-band power (z-score)
%cfg_def.methods = 'raw'; % not used?
cfg_def.merge_thr = 0.02; % merge ivs closer together than this (s)
cfg_def.minlen = 0.02; % iv minimum length (s)

cfg_def.viewMUA = 0;
cfg_def.useMUA = 0;
cfg_def.overlap_thr = 0.01; %in seconds
cfg_def.display_filter = 0;

cfg = ProcessConfig2(cfg_def,cfg_in);

%% Get SWR from lfp power
cfg.methods = 'raw';
% bandpass lfp
l_filter = FilterLFP(cfg,lfp);

% obtain power and z-score
l_power = LFPpower([],l_filter);
l_zscore = zscore_tsd(l_power);

% threshold events
spwr_iv = TSDtoIV(cfg,l_zscore);

% to each event, add a field with the max z-scored power
cfg_temp = [];
cfg_temp.method = 'max'; % 'min', 'mean'
cfg_temp.label = 'maxz'; % what to call this in iv, i.e. usr.label
spwr_iv = AddTSDtoIV(cfg_temp,spwr_iv,l_zscore);

% select only those events of >5 z-scored power
cfg_temp.threshold = cfg.max_thr;
spwr_iv = SelectIV(cfg_temp,spwr_iv);
spwr_iv_out = spwr_iv;

%% get multi-unit activity

% obtain multi-unit activity and z-score
cfg.tvec = lfp.tvec;
mua = getMUA(cfg,S); 
mua_z = zscore_tsd(mua);

% threshold mua
muaIV = TSDtoIV(cfg,mua_z);

% to each event, add a field with the max z-scored power
cfg_temp = [];
cfg_temp.method = 'max'; % 'min', 'mean'
cfg_temp.label = 'MUA_maxz'; % what to call this in iv, i.e. usr.label
muaIV = AddTSDtoIV(cfg_temp,muaIV,mua_z);

% select only those events of >5 z-scored power
cfg_temp = [];
cfg_temp.threshold = cfg.max_thr;
muaIV = SelectIV(cfg_temp,muaIV);

if cfg.viewMUA
    PlotTSDfromIV([],muaIV,mua_z);
end

%% Combine lfp power and MUA
% Keep those time intervals in where MUA is over a threshold and where sharp-wave ripples
% occur within some timewindow

if cfg.useMUA
    numIV = length(muaIV.tstart);
    keep_idx = [];

    for iIV = 1:numIV
        curr_tstart = muaIV.tstart(iIV) - cfg.overlap_thr;
        curr_tend = muaIV.tend(iIV) + cfg.overlap_thr;

        check1 = any(find(curr_tstart < spwr_iv.tstart & spwr_iv.tstart < curr_tend));
        check2 = any(find(curr_tstart < spwr_iv.tend & spwr_iv.tend < curr_tend));

        if check1 || check2
            keep_idx = cat(1,keep_idx,iIV);
        end

    end

    keep_tstart = muaIV.tstart(keep_idx);
    keep_tend = muaIV.tend(keep_idx);

    spwr_iv_out = iv(keep_tstart,keep_tend);
end





end