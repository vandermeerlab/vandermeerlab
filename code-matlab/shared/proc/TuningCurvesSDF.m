function tc = TuningCurvesSDF(cfg_in,sdf,tuning_var)
% function tc = TuningCurvesSDF(cfg_in,sdf,tuning_var)
%
% computes tuning curves for spike density functions in sdf relative to
% variable(s) tuning_var (really just average values for binned samples)
%
% INPUTS:
%
% sdf: tsd with spike density functions [nCells x nSamples]
% tuning_var: tsd with [nDim x nSamples] tuning variable (e.g. position)
%
% OUTPUTS:
%
% tc.tc: nCells x nBins tuning curves (average firing rate)
% tc.occ: 1 x nBins occupancy (sample count)
%
% CFG:
%
% cfg.binEdges: {1 x nDim} cell array with bin edges (default obtained from
%  min, max in data)
% cfg.occ_dt = 1/30; % time corresponding to each occupancy sample
%
% MvdM 2016-05-07

cfg_def = [];
cfg_def.binEdges = linspace(min(tuning_var.data(1,:)),max(tuning_var.data(1,:)),101); % 100 bins s
cfg_def.occ_dt = 1/30;

cfg = ProcessConfig(cfg_def,cfg_in);

%
if ~CheckTSD(sdf) | ~CheckTSD(tuning_var)
   error('Incorrect TSD.'); 
end

if size(sdf.data,2) ~= size(tuning_var.data,2)
   error('Unequal numbers of samples.'); 
end

if size(tuning_var.data,1) > 1
   error('>1D not yet implemented.'); 
end

%
nCells = size(sdf.data,1);
nBins = length(cfg.binEdges)-1;

tc.tc = nan(nCells,nBins);
tc.occ = nan(1,nBins);

for iC = 1:nCells
    
    [~,bin] = histc(tuning_var.data(1,:),cfg.binEdges);
    bin(bin == nBins + 1) = nBins; % data points that are on last edge go with previous bin
    
    for iB = 1:nBins
    
        this_bin_data = sdf.data(iC,bin == iB);
        tc.tc(iC,iB) = nanmean(this_bin_data);
        tc.occ(iB) = length(this_bin_data); % only need to do this once, but oh well
        
    end
    
end

tc.occ = tc.occ .* cfg.occ_dt;