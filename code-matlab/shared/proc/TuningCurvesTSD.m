function [tc, nSamples] = TuningCurvesTSD(cfg_in, target_tsd, tuning_tsd)
% function [tc, nSamples] = TuningCurvesTSD(cfg_in, target_tsd, tuning_tsd)
%
% computes tuning curves (binned averages) for dependent variable target_tsd relative to
% variable(s) tuning_tsd
%
% examples: average theta power as a function of position, dopamine level as a
% function of speed
%
% INPUTS:
%
% target_tsd: variable of interest; the data to be binned , [1 x nSamples] tsd
% tuning_var: variable used for binning (e.g. position), [nDim x nSamples] tsd
%
% OUTPUTS:
%
% tc: [nBins x 1] (1D) or [nBinsX x nBinsY] (2D) matrix of average target_tsd values binned by tuning_tsd 
% nSamples: target_tsd counts for each tc bin
%
% CFG:
%
% cfg.binEdges = []; % {1 x nDim} cell array with bin edges to apply to tuning_tsd (default obtained from
%  min, max in tuning_tsd.data with cfg.nBins bins)
% cfg.nBins = [100 100]; % only used if cfg.binEdges is empty, use this many bins
%  between min and max for each dimension (can be 1D or 2D) 
% cfg.minOcc = 1; % minimum occupancy (in samples)
%
%%% EXAMPLE: compute average Y position for X, Y bins
%
% pos = LoadPos([]);
%
% posX = pos; posX.data = posX.data(1,:); posX.label = posX.label(1);
% posY = pos; posY.data = posY.data(2,:); posY.label = posY.label(2);
%
% [tc, nSamples] = TuningCurvesTSD([], cfg, posY, pos);
% imagesc(tc);
%
% MvdM 2022

% preamble

if ~CheckTSD(target_tsd)
   error('target_tsd input is not a well-formed tsd.'); 
end

if ~CheckTSD(tuning_tsd)
   error('tuning_tsd input is not a well-formed tsd.'); 
end

cfg_def = [];
cfg_def.binEdges = [];
cfg_def.nDim = size(tuning_tsd.data, 1); 

if ~isfield(cfg_in, 'binEdges') % set up default bins
    
    cfg_def.nBins = [100 100];
    if isfield(cfg_in, 'nBins')
        cfg_def.nBins = cfg_in.nBins;
    end
    
    for iDim = 1:cfg_def.nDim
        
        mn = min(tuning_tsd.data(iDim, :));
        mx = max(tuning_tsd.data(iDim, :));
        
        cfg_def.binEdges{iDim} = linspace(mn, mx, cfg_def.nBins(iDim) + 1);
    end
    
else % bin edges defined, need to discover nBins
    
    for iDim = 1:cfg_def.nDim
        cfg_def.nBins(iDim) = length(cfg_in.binEdges{iDim}) - 1;
    end
    
end

cfg_def.minOcc = 1;

cfg = ProcessConfig(cfg_def, cfg_in);

% work
switch cfg.nDim
    case 1
      
        [~, bin] = histc(tuning_tsd.data, cfg.binEdges{1}); % find out bin idx for tuning data
        bin(bin == cfg.nBins(1) + 1) = cfg.nBins(1); % data points that are on last edge go with previous bin
        
        for iB = cfg.nBins(1):-1:1
            
            this_bin_data = target_tsd.data(bin == iB);
            
            tc(iB) = nanmean(this_bin_data);
            nSamples(iB) = length(this_bin_data);
            if nSamples(iB) < cfg.minOcc, tc(iB) = NaN; end 
            
        end
        
    case 2
        
        [~, ~, ~, bin1, bin2] = histcounts2(tuning_tsd.data(1,:), tuning_tsd.data(2,:), ...
            cfg.binEdges{1}, cfg.binEdges{2});
        bin1(bin1 == cfg.nBins(1) + 1) = cfg.nBins(1); bin2(bin2 == cfg.nBins(2) + 1) = cfg.nBins(2);
        
        for iB1 = cfg.nBins(1):-1:1 % bin1 are idxs of those samples in ea
          
            for iB2 = cfg.nBins(2):-1:1
                
                this_bin_data = target_tsd.data(bin1 == iB1 & bin2 == iB2);
                tc(iB1, iB2) = nanmean(this_bin_data);
                nSamples(iB1, iB2) = length(this_bin_data);
                
            end
                        
        end
        
    otherwise
        error('More than 2 tuning dimensions is not yet implemented.');
end




