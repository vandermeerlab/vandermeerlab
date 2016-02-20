function Q = MakeQfromS2(cfg_in,S,iv_in)
%MAKEQFROMS2 Obtain spike counts for S in time bins defined by iv.
%
%function Q = MAKEQFROMS2(cfg,S,iv) makes a "Q-matrix" (spike counts in 
% [nCells x nTimeBins] matrix) from ts with spike data. nCells =
% length(S.t) and nTimeBins = length(iv.tstart).
%
%  A Q-matrix is organized such at each row represents a neural unit and
%  each column represents a timebin that units are potentially active in.
%  Thus, each element of the Q-matrix is a count of the number of spikes 
%  emitted by the correpsonding neural unit in the corresponding time bin. 
%
%    INPUTS:
%       S: ts datatype containing spike times
%       iv_in: iv datatype containing intervals that will be the basis for
%              binning spikes.
%
%    OUTPUTS:
%
%       Q: tsd with fields
%          .data = [nCells x nTimeBins] Q-matrix containing spike counts
%          .tvec timebin centers
%          .label {1 x nCells} identifies which row of the Q-matrix
%               corresponds to which neural unit from S.
%
%    CONFIG OPTIONS:
%       cfg.win = []; Redefine a new window around the center of the
%              intervals. Ex: cfg.win = 0.1 is asking that the intervals be
%              modified so that they have a duration of 100 ms.
%       cfg.smooth = []; Smoothing window shape. If empty [], no smoothing
%              is done; otherwise shape options are 'gauss' and 'boxcar'.
%       cfg.gausswin_size = 1; % in seconds; only used if cfg.smooth is
%              'gauss'.
%       cfg.gausswin_sd = 0.02; % in seconds, applies only if cfg.smooth is
%              'gauss.
%       cfg.boxcar_size = 5; % in nBins, applies only if cfg.smooth is
%              'boxcar'.
%
% MvdM 2014-08-21
% aacarey Dec 2015, version 2

cfg_def = [];
cfg_def.win = [];
cfg_def.smooth = []; % [], 'gauss', 'boxcar'
cfg_def.gausswin_size = 1; % in seconds; only used if cfg.smooth not empty
cfg_def.gausswin_sd = 0.02; % in seconds
cfg_def.boxcar_size = 5; % in nBins
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% check inputs
if isempty(S.t)
    Q = tsd([],[]); warning('S.t is empty!');
    return;
end

nCells = length(S.t);
nTimeBins = length(iv_in.tstart);

% preallocate space
Qdata = nan(nCells,nTimeBins);

if ~isempty(cfg.win); ivCenters = IVcenters(iv_in); end % compute once
for iIV = length(iv_in.tstart):-1:1
    if isempty(cfg.win) % use exact boundaries
        cfg.tvec_edges = [iv_in.tstart(iIV) iv_in.tend(iIV)];
        
    else % redfine boundaries according to cfg.win
        center_t = ivCenters(iIV);
        cfg.tvec_edges = [center_t-cfg.win/2 center_t+cfg.win/2];
    end
    
    cfg.dt = diff(cfg.tvec_edges);
    cfg.verbose = 0;
    qtemp = MakeQColumn(cfg,S); % helper function, from original MakeQfromS
    tvec(iIV) = qtemp.tvec;
    Qdata(:,iIV) = qtemp.data;
    
end % of events

% make output
Q = tsd(tvec,Qdata,qtemp.label);
Q.usr = S.usr;

% housekeeping
Q.cfg.history.mfun = cat(1,Q.cfg.history.mfun,mfun);
Q.cfg.history.cfg = cat(1,Q.cfg.history.cfg,{cfg});

%~~~~~~~~~ helper hack ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function Qcolumn = MakeQColumn(cfg,S)
        
        % if smoothing, set up kernel
        if ~isempty(cfg.smooth)
            switch cfg.smooth
                case 'gauss'
                    gauss_window = cfg.gausswin_size./cfg.dt; % 1 second window
                    gauss_SD = cfg.gausswin_sd./cfg.dt; % 0.02 seconds (20ms) SD
                    gk = gausskernel(gauss_window,gauss_SD)'; %gk = gk./cfg.dt; % normalize by binsize
                case 'boxcar'
                    gk = ones(1,cfg.boxcar_size);
            end
        end
        
        % construct Q-matrix
        for iC = length(S.t):-1:1
            
            spk_t = S.t{iC};
            Qcolumn(iC,:) = trim_histc(histc(spk_t,cfg.tvec_edges));
            
            if ~isempty(cfg.smooth)
                Qcolumn(iC,:) = conv2(Qcolumn(iC,:),gk,'same'); % convolve with smoothing kernel
            end
            label(iC) = S.label(iC);
        end
        
        Qcolumn = tsd(mean(cfg.tvec_edges),Qcolumn,label);
        
    end % of MakeQfromS1
end % of MakeQfromS2
