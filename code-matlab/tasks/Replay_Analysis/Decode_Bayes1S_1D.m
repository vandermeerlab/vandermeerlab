function iv_out = Decode_Bayes1S_1D(cfg_in,S,tc)
% function iv_out = Decode_Bayes1S_1D(cfg_in,S,tc)
%
% one-step (naive) Bayesian decoder for one-dimensional variables
%
% INPUTS:
%
% S: ts with spike data (e.g. from LoadSpikes())
% tc: tuning curves (e.g. from MakeTC())
%
% *cfg_in must contain the evt field with time intervals to decode
%
% OUTPUTS:
%
% iv_out: intervals with various usr fields of interest
% .Decode_Bayes1S_1D: posterior (decoded) probability p(variable=value | spikes)
% .Decode_Bayes1S_1D_Qmat: Q-matrix
% .Decode_Bayes1S_1D_pshuf: posterior for each of nShuffles shuffles
%
% OPTIONS:
%
% cfg_def.dt = 0.025; % time step (in seconds)
% cfg_def.smooth = 0; % % [], 'gauss', 'boxcar'
% cfg_def.gausswin_size = 1; % Gaussian smoothing kernel for Q-matrix, width in seconds
% cfg_def.gausswin_sd = 0.02; % smoothing kernel SD (in seconds)
% cfg_def.boxcar_size = 5; % boxcar size in bins
% cfg_def.trim = 1; % trim leading and trailing zeros in Q-matrix (for speed later)
% cfg_def.returnQ = 0; % return full Q-matrix in output
% cfg_def.regularizeTC = 0; % regularize tuning curves by replacing zero values with this
% cfg_def.shuffle_mode = 0; % if ~= 0, perform shuffles (see help ShuffleTS to see the options)
% cfg_def.nShuffles = 100; % number of shuffles to perform
%
% MvdM Feb 2015

cfg_def = [];
cfg_def.dt = 0.025; % time step (in seconds)
cfg_def.smooth = 0; % smooth rows of Q-matrix with Gaussian kernel if set to 1
cfg_def.gausswin_size = 1; % Gaussian smoothing kernel for Q-matrix, width in seconds
cfg_def.gausswin_sd = 0.02; % smoothing kernel SD (in seconds)
cfg_def.boxcar_size = 5; % boxcar size in bins
cfg_def.trim = 1; % trim leading and trailing zeros in Q-matrix (for speed later)
cfg_def.returnQ = 0; % return full Q-matrix in output
cfg_def.regularizeTC = 0; % regularize tuning curves by replacing zero values with this
cfg_def.shuffle_mode = 0; % if ~= 0, perform shuffles (see help ShuffleTS to see the options)
cfg_def.nShuffles = 100; % number of shuffles to perform

cfg = ProcessConfig2(cfg_def,cfg_in);
mfun = mfilename;

% some checks on input

if ~isfield(cfg,'evt')
    error('No cfg.evt field specified.');
end

%
if size(tc,2) ~= length(S.t)
   error('size(2) of tc (%d) is not consistent with number of cells in S.t (%d)',size(tc,2),length(S.t)); 
end

iv_out = cfg.evt;

%
nUsr = 1;
if isfield(iv_out,'usr')
    nUsr = length(iv_out.usr) + 1;
end

% regularize TCs
if cfg.regularizeTC ~= 0
   tc(tc == 0) = cfg.regularizeTC; 
end

% main loop

nIV = length(cfg.evt.tstart);
for iIV = nIV:-1:1 % NOTE: just decoding the whole thing and restricting later may be faster...

    fprintf('Event %d/%d...\n',iIV,nIV);
    
    % make Q-matrix (spike counts as nCells x nTimeBins) for this iv
    
    cfg_Q = cfg;
    
    % define time base for Q-matrix
    cfg_Q.tvec_edges = cfg.evt.tstart(iIV):cfg.dt:cfg.evt.tend(iIV);
    
    % call dec_function here to create "actually observed" decoding
    
    [p,Q] = dec_function(cfg_Q,S,tc);
    
    iv_out.usr(nUsr).data{iIV} = p;
    
    if cfg.returnQ
        iv_out.usr(nUsr+1).data{iIV} = Q;
    end
    
    % if shuffles specified, loop here to create shuffled S's, call
    % dec_function again
    
    if cfg.shuffle_mode ~= 0
        
        cfg_shuf = []; cfg_shuf.mode = cfg.shuffle_mode;
        %cfg_shuf.t0 = cfg.evt.tstart(iIV); cfg_shuf.t1 = cfg.evt.tend(iIV); 
        % start and end times of shuffle should match edges used to make
        % Q and p matrices -- this is very important!
        tvec_centers = iv_out.usr(nUsr).data{iIV}.tvec;
        
        if isempty(tvec_centers)
            iv_out.usr(nUsr+2).data{iIV} = [];
            continue;
        end
        
        cfg_shuf.t0 = tvec_centers(1)-cfg_Q.dt/2;
        cfg_shuf.t1 = tvec_centers(end)+cfg_Q.dt/2;
        
        cfg_Qshuf = cfg_Q; cfg_Qshuf.trim = 0; % don't trim shuffled spikes!
        
        for iShuf = cfg.nShuffles:-1:1
                   
            S_shuf = ShuffleTS(cfg_shuf,S);
            p = dec_function(cfg_Qshuf,S_shuf,tc);
            
            iv_out.usr(nUsr+2).data{iIV}(iShuf,:,:) = p;
            
        end
    end

    
end % of loop over ivs

iv_out.usr(nUsr).label = 'Decode_Bayes1S_1D';

if cfg.returnQ
    iv_out.usr(nUsr+1).label = 'Decode_Bayes1S_1D_Qmat';
end

if cfg.shuffle_mode ~= 0
    iv_out.usr(nUsr+2).label = 'Decode_Bayes1S_1D_pshuf';
end

% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_out.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_out.cfg.history.cfg,{cfg});


function [p,Q] = dec_function(cfg_Q,S,tc)
% this is the actual decoding function

Q = MakeQfromS(cfg_Q,S);

nActiveNeurons = sum(Q.data > cfg_Q.regularizeTC); % this affects the regression quite a bit...

% decode
tvec_centers = cfg_Q.tvec_edges(1:end-1)+cfg_Q.dt/2;

if cfg_Q.trim
    keep_idx = find(nActiveNeurons,1,'first'):find(nActiveNeurons,1,'last');
    nActiveNeurons = nActiveNeurons(keep_idx);
    tvec_centers = tvec_centers(keep_idx);
    Q.data = Q.data(:,keep_idx);
    Q.tvec = Q.tvec(keep_idx);
end

if isempty(tvec_centers)
    p.data = []; p.tvec = []; Q.data = []; Q.tvec = [];
    return;
end;

nBins = size(tc,1);
len = length(tvec_centers);
p = nan(length(tvec_centers),nBins);

occUniform = ones(nBins,1);

for iB = 1:nBins
    
    tempProd = nansum(log(repmat(tc(iB,:)',1,len).^Q.data));
    tempSum = exp(-cfg_Q.dt*nansum(tc(iB,:),2));
    
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
    
end

p = p./repmat(sum(p,2),1,nBins); % normalize
p(nActiveNeurons < 1,:) = 0;

p = tsd(tvec_centers,p);