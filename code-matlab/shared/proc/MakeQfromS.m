function Q = MakeQfromS(cfg_in,S)
% function Q = MakeQfromS(cfg,S)
%
% Make "Q-matrix" (spike counts in [nCells x nTimeBins] matrix) from ts with spike data
%
% INPUTS:
%
% S: ts with spikes
%
% CONFIG OPTIONS:
%
% cfg_def = [];
% cfg_def.dt = 0.05;
% cfg_def.smooth = 0;
% cfg_def.gausswin_size = 1; % in seconds; only used if cfg.smooth = 1
% cfg_def.gausswin_sd = 0.02; % in seconds
% cfg_def.tvec_edges = [];
% 
% OUTPUTS:
%
% Q: [nCells x nTimeBins] tsd with spike counts
%
% MvdM 2014-08-21

cfg_def = [];
cfg_def.dt = 0.05;
cfg_def.smooth = []; % [], 'gauss', 'boxcar'
cfg_def.gausswin_size = 1; % in seconds; only used if cfg.smooth = 1
cfg_def.gausswin_sd = 0.02; % in seconds
cfg_def.boxcar_size = 5; % in nBins
cfg_def.tvec_edges = [];

cfg = ProcessConfig2(cfg_def,cfg_in);
mfun = mfilename;

% check inputs
if isempty(S.t)
    Q = tsd([],[]); warning('S.t is empty!');
    return;
end

% assemble tvecs
if isempty(cfg.tvec_edges)
    cfg.tvec_edges = firstSpike(S):cfg.dt:lastSpike(S);
else
    cfg.dt = median(diff(cfg.tvec_edges));
end

tvec_centers = cfg.tvec_edges(1:end-1)+cfg.dt/2;

% if smoothing, set up kernel
if ~isempty(cfg.smooth)
    switch cfg.smooth
        case 'gauss'
            gauss_window = cfg_def.gausswin_size./cfg.dt; % 1 second window
            gauss_SD = cfg_def.gausswin_sd./cfg.dt; % 0.02 seconds (20ms) SD
            gk = gausskernel(gauss_window,gauss_SD)'; %gk = gk./cfg.dt; % normalize by binsize
        case 'boxcar'
            gk = ones(1,cfg.boxcar_size);
    end
end

% construct Q-matrix
for iC = length(S.t):-1:1
    
    spk_t = S.t{iC};
    Q(iC,:) = trim_histc(histc(spk_t,cfg.tvec_edges));
    
    if ~isempty(cfg.smooth)
        Q(iC,:) = conv2(Q(iC,:),gk,'same'); % convolve with smoothing kernel
    end
    
end

Q = tsd(tvec_centers,Q);

% housekeeping
Q.cfg.history.mfun = cat(1,Q.cfg.history.mfun,mfun);
Q.cfg.history.cfg = cat(1,Q.cfg.history.cfg,{cfg});
