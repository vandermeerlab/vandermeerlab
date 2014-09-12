function p_tsd = DecodeZ(cfg_in,Q,tc)
% function p = DecodeZ(cfg,Q,tc)
%
% perform Bayesian decoding on input Q-matrix using encoding model
% specified by tuning curves in tc
%
%
% INPUTS:
%
% Q: tsd with [nCells x nTimeBins double] Q-matrix of spike counts (firing rates), from MakeQfromS()
% tc: [nCells x nVarBins double] tuning curves, from TuningCurves()
%
% OUTPUTS:
%
% p: tsd with [nTimeBins x nVarBins double] matrix with posterior p(x|s)
%
% CFG OPTIONS:
%
% cfg.noSpikesInBin = 'zeros'; % what to do with time bins without spikes
% cfg.mode = 'oneStep'; % 'twoStep'; algorithm to use, see Zhang et al. (1998)
% cfg.kernel = gausskernel(10,1); % kernel to use for model-based
%   prior (only used if cfg.mode = 'twoStep')
% cfg.resetTimes = []; % times in which model-based prediction is
%   discarded so that only likelihood is used (only used if cfg.mode =
%   'twoStep')
%
% MvdM 2014-08-22 initial version

cfg = [];
cfg.noSpikesInBin = 'zeros';
cfg.kernel = gausskernel(10,1);
cfg.resetTimes = [];
cfg.mode = 'oneStep';

ProcessConfig;
mfun = mfilename;

%
len = length(Q.tvec);
nBins = size(tc,2);
binsize = median(diff(Q.tvec)); % NOTE: should have a stronger check for this

occUniform = ones(1,nBins)./nBins;
p = nan(len,nBins);

% main decoding loop
for iB = 1:nBins
    tempProd = nansum(log(repmat(tc(:,iB),1,len).^Q.data));
    tempSum = exp(-binsize*nansum(tc(:,iB)));
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
end
p = p./repmat(sum(p,2),1,nBins); % renormalize

% deal with bins without any spikes
nActiveNeurons = sum(Q.data > 0);

switch cfg.noSpikesInBin
    case 'zeros'
        p(nActiveNeurons < 1,:) = 0;
    case 'nans'
        p(nActiveNeurons < 1,:) = NaN;
end

%
p_tsd = tsd(Q.tvec,p);
p_tsd.usr.nActiveNeurons = nActiveNeurons;

% housekeeping
p_tsd.cfg.history.mfun = cat(1,p_tsd.cfg.history.mfun,mfun);
p_tsd.cfg.history.cfg = cat(1,p_tsd.cfg.history.cfg,{cfg});