function [x,y,Timestamps] = xy_targets_chunks(cfg_in)
% function [x,y,Timestamps] = xy_targets_chunks(cfg_in)
%
% wrapper function for Neuralynx ExtractFromTargets.m function
%
% assumes VT1.nvt available file in working directory
%
% chunks target extraction into specified size for speedup
%
% OUTPUTS:
%
% x: nSamples x nTargets x-coordinates (pixels)
% y: nSamples x nTargets y-coordinates (pixels)
% Timestamps: RAW Neuralynx timestamps
%
% INPUTS:
%
% cfg_default.targets = []; % targets to extract, e.g. [] (all), 1, 1:5
% cfg_default.mode = 'raw'; % 'avg','raw' % avg averages over multiple targets
% cfg_default.chunksize = 5000; % chunk size in samples
% cfg_default.var_tolerance = 500; % do not average if variance exceeds this (pixels)
%
%
% MvdM 2015-04-26 based on EI initial version

cfg_default.targets = []; % targets to extract, e.g. [] (all), 1, 1:5
cfg_default.mode = 'raw'; % 'avg','raw'
cfg_default.chunksize = 5000;
cfg_default.var_tolerance = 500;

cfg = ProcessConfig2(cfg_default,cfg_in);

% load raw data
[Timestamps, X, Y, Angles, Targets, Points, Header] = ...
    Nlx2MatVT(FindFile('*VT1.nvt'), [1 1 1 1 1 1], 1, 1, []);

% report some stuff about outputs
TargetRowAllZeros = all(Targets' == 0);
nTargets = sum(~TargetRowAllZeros);
fprintf('%d targets found in full file (non-zero rows in Targets).\n',nTargets);

% set up 
nBins = length(Targets);
chunks = [1:cfg.chunksize-1:nBins, nBins];

% initialize output variables x, y (depends on number of targets!)
switch cfg.mode
    
    case 'raw'
        if isempty(cfg.targets) % use all
            NoutTargets = nTargets;
        else
            NoutTargets = length(cfg.targets);
        end
        
    case 'avg'
        NoutTargets = 1;
end
x = nan(nBins,NoutTargets);
y = nan(nBins,NoutTargets);

% loop for target extraction
for iC = 1:length(chunks)-1
    fprintf('Processing chunk %d of %d...\n',iC,length(chunks));
    small_target_matrix = Targets(:,chunks(iC):chunks(iC+1));
    [x_temp, y_temp, color_temp, valid_targets_temp] = ...
        ExtractFromTargets(small_target_matrix);
    x_temp(x_temp == 0) = nan;
    y_temp(y_temp == 0) = nan;
        
    %
    nTargetsFound = max(valid_targets_temp);
    if nTargetsFound ~= size(x_temp,2)
        error('Inconsistent number of targets found.');
    end
    
    fprintf('%d targets found in chunk.\n',nTargetsFound);
    
    % select targets
    if ~isempty(cfg.targets)
    
        if max(cfg.targets) > nTargetsFound
            
            warning('Specified target not found.');
            
            x_temp = nan(size(x_temp));
            y_temp = nan(size(x_temp));
            
        else
            
            x_temp = x_temp(:,cfg.targets);
            y_temp = y_temp(:,cfg.targets);
            
        end
    
    end
    
    %
    switch cfg.mode
        
        case 'raw'
            
            x(chunks(iC):chunks(iC)+length(x_temp(:,1))-1,1:size(x_temp,2)) = x_temp;
            y(chunks(iC):chunks(iC)+length(y_temp(:,1))-1,1:size(y_temp,2)) = y_temp;
        
        case 'avg'
            
            % Averaging over targets when varience isn't too high, otherwise nan.
            x_avg = nanmean(x_temp');
            x_var = nanvar(x_temp');
            x_avg(x_var > cfg.var_tolerance) = nan;
            x(chunks(iC):chunks(iC)+length(x_temp(:,1))-1,:) = x_avg;
            
            y_avg = nanmean(y_temp');
            y_var = nanvar(y_temp');
            y_avg(y_var > cfg.var_tolerance) = nan;
            y(chunks(iC):chunks(iC)+length(y_temp(:,1))-1,:) = y_avg;
            
    end
    
%     % Using the first target only.
%     x(chunks(i):chunks(i)+length(x_temp(:,1))-1,:) = x_temp(:,1);
%     y(chunks(i):chunks(i)+length(y_temp(:,1))-1,:) = y_temp(:,1);

end


