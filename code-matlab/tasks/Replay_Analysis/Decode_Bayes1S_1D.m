function iv_out = Decode_Bayes1S_1D(cfg_in,S,tc)
% function iv_out = Decode_Bayes1S_1D(cfg_in,S,tc)
%


cfg_def = [];
cfg_def.dt = 0.005; % in ms

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

% main loop

nIV = length(cfg.evt.tstart);
for iIV = nIV:-1:1 % NOTE: just decoding the whole thing and restricting later may be faster...

    % make Q-matrix (spike counts as nCells x nTimeBins) for this iv
    
    cfg_Q = [];
    
    % define time base for Q-matrix
    cfg_Q.tvec_edges = cfg.evt.tstart(iIV):cfg.dt:cfg.evt.tend(iIV);
    
    Q = MakeQfromS(cfg_Q,S);
    
    nActiveNeurons = sum(Q.data > 0);
    
    % optionally, could smooth here (should make smoothQMat() function...)
    
    % decode
    tvec_centers = cfg_Q.tvec_edges(1:end-1)+cfg.dt/2;
    
    len = length(tvec_centers);
    
    nBins = size(tc,1);
    p = nan(length(tvec_centers),nBins);
    
    occUniform = ones(nBins,1);
    
    for iB = 1:nBins
        
        tempProd = nansum(log(repmat(tc(iB,:)',1,len).^Q.data));
        tempSum = exp(-cfg.dt*nansum(tc(iB,:),2));
        
        p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
        
    end
    
    p = p./repmat(sum(p,2),1,nBins); % normalize
    p(nActiveNeurons < 1,:) = 0;
    
    p = tsd(tvec_centers,p);
    
    iv_out.usr(nUsr).data{iIV} = p;
    
end % of loop over ivs

iv_out.usr(nUsr).label = 'Decode_Bayes1S_1D';



% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_out.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_out.cfg.history.cfg,{cfg});
