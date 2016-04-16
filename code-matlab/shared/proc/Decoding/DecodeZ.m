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
% cfg.noSpikesInBin = 'zeros'; % {'zeros','nans'}, what to put in time bins
%   below threshold
% cfg_def.nMinSpikes = 1; % minimum number of spike (sum) for bin to be
%   included
%
% NOTE: assumes Q.data is in format [nCells x nTimeBins]
%
% MvdM 2014-08-22 initial version, later simplified to include only
% one-step

cfg_def = [];
cfg_def.noSpikesInBin = 'nans';
cfg_def.nMinSpikes = 1;

cfg = ProcessConfig(cfg_def,cfg_in);
mfun = mfilename;

%
len = length(Q.tvec);
nBins = size(tc,2);
binsize = median(diff(Q.tvec)); % NOTE: should have a stronger check for this

occUniform = ones(1,nBins)./nBins;
p = nan(len,nBins);

% main decoding loop -- note, use a version of nansum that returns nan when summing over
% nans only!
for iB = 1:nBins
    tempProd = nansum(log(repmat(tc(:,iB),1,len).^Q.data));
    tempSum = exp(-binsize*nansum(tc(:,iB)));
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
end
p = p./repmat(nansum(p,2),1,nBins); % renormalize

% deal with bins without any spikes
nActiveNeurons = sum(Q.data > cfg.nMinSpikes);

switch cfg.noSpikesInBin
    case 'zeros'
        p(nActiveNeurons < 1,:) = 0;
    case 'nans'
        p(nActiveNeurons < 1,:) = NaN;
end

%
p_tsd = tsd(Q.tvec,p');
p_tsd.usr.nActiveNeurons = nActiveNeurons;

if ~CheckTSD(p_tsd)
   error('DecodeZ: malformed TSD.'); 
end

% housekeeping
p_tsd.cfg.history.mfun = cat(1,p_tsd.cfg.history.mfun,mfun);
p_tsd.cfg.history.cfg = cat(1,p_tsd.cfg.history.cfg,{cfg});