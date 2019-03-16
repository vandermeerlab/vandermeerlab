function p_tsd = DecodeZ2(cfg_in,Q,tc)
% function p = DecodeZ2(cfg,Q,tc)
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
% cfg_def.excludeMethod = 'frate'; % 'frate', 'nNeurons'; criterion for excluding time bin
% cfg_def.nMinNeurons = 1; % if cfg.excludeMethod = 'nNeurons', specify how many neurons need to be active
% cfg_def.nMinSpikes = 1; % threshold for counting neuron as active (if cfg.excludeMethod = 'nNeurons') or bin as active (if cfg_def.excludeMethod = 'frate')
%
% NOTE: assumes Q.data is in format [nCells x nTimeBins]
%
% MvdM 2014-08-22 initial version, later simplified to include only
% one-step
% Elyot 14 Mar 2019 - speed improvement, removed loops and repmats and
% vectorized everything. Had to do some wrangling to preserve existing NaN
% behaviour. May differ in output by a few parts per million (usually much
% less) because of numerical instability... basically
% exp(sum(log(a^b))) has been replaced with exp(sum(b*log(a))) and the
% whole computation can be a bit unstable with values extremely close to 0

cfg_def = [];
cfg_def.noSpikesInBin = 'nans';
cfg_def.excludeMethod = 'frate'; % 'frate', 'nNeurons'; criterion for excluding time bin
cfg_def.nMinNeurons = 1; % if cfg.excludeMethod = 'nNeurons', specify how many neurons need to be active
cfg_def.nMinSpikes = 1; % threshold for counting neuron as active (if cfg.excludeMethod = 'nNeurons') or bin as active (if cfg_def.excludeMethod = 'frate')

cfg = ProcessConfig(cfg_def,cfg_in);
mfun = mfilename;

%
len = length(Q.tvec);
nBins = size(tc,2);
binsize = median(diff(Q.tvec)); % NOTE: should have a stronger check for this

occUniform = ones(1,nBins)./nBins;

% main decoding loop -- note, use a version of nansum that returns nan when summing over
% nans only!
tempQdata = Q.data;
nans = isnan(tempQdata);
tempQdata(nans)=0;

templogtc =log(tc);
tcnans = isnan(templogtc);
templogtc(tcnans) = 0;
tcinfs = isinf(templogtc);
templogtc(tcinfs) = -realmax; % Must not be -inf because we need it to not give nan when multiplied by zero.
tempProd = tempQdata' * templogtc;

tempProd(all(nans),:) = NaN;

%tempSum = exp(-binsize*sum(tc,1));
tempSum = exp(-binsize*nansum(tc,1));
p = exp(tempProd).*tempSum.*occUniform;


%for iB = 1:nBins
%    %tempProd = nansum(log(repmat(tc(:,iB),1,len).^Q.data));
%    
%    tempSum = exp(-binsize*nansum(tc(:,iB)));
%    p(:,iB) = exp(tempProd(iB,:))*tempSum*occUniform(iB);
%end

%p = p./repmat(nansum(p,2),1,nBins); % renormalize

p = p./(nansum(p,2));
%p = p./(sum(p,2));

% deal with bins without any spikes
switch cfg.excludeMethod
case 'nNeurons' % number of neurons with at least this activity level
    nActiveNeurons = sum(Q.data >= cfg.nMinSpikes);
    toss_idx = nActiveNeurons < cfg.nMinNeurons;
case 'frate'
    nActiveNeurons = sum(Q.data); % this is an overall activity level; not how many individual neurons pass an activity threshold
    toss_idx = nActiveNeurons < cfg.nMinSpikes;
end

switch cfg.noSpikesInBin
case 'zeros'
    p(toss_idx,:) = 0;
case 'nans'
    p(toss_idx,:) = NaN;
end

%
p_tsd = tsd(Q.tvec,p');
p_tsd.usr.nActiveNeurons = nActiveNeurons;
p_tsd.usr.nActiveNeuronsPassed = (nActiveNeurons >= cfg.nMinSpikes);

if ~CheckTSD(p_tsd)
error('DecodeZ: malformed TSD.'); 
end

% housekeeping
p_tsd.cfg.history.mfun = cat(1,p_tsd.cfg.history.mfun,mfun);
p_tsd.cfg.history.cfg = cat(1,p_tsd.cfg.history.cfg,{cfg});