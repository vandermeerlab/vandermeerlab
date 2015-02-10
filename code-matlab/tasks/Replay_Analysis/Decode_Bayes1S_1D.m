function iv_out = Decode_Bayes1S_1D(cfg_in,S,tc)
% function iv_out = Decode_Bayes1S_1D(cfg_in,S,tc)
%


cfg_def = [];
cfg_def.dt = 0.025; % in ms
cfg_def.smooth = 1;
cfg_def.gwin = 1; % s
cfg_def.gsd = 0.02; % s
cfg_def.returnQ = 0;
cfg_def.regularizeTC = 0;
cfg_def.shuffle_mode = 0;
cfg_def.nShuffles = 100;

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
        cfg_shuf.t0 = cfg.evt.tstart(iIV); cfg_shuf.t1 = cfg.evt.tend(iIV); 
        
        for iShuf = cfg.nShuffles:-1:1
                   
            S_shuf = ShuffleTS(cfg_shuf,S);
            p = dec_function(cfg_Q,S_shuf,tc);
            
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
len = length(tvec_centers);

nBins = size(tc,1);
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