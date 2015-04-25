function [ylims] = PlotSpikeRaster2(cfg_in,S_in)
%% Plot spikes
% Converting the data to (x,y) coordinates is significantly faster to plot than calling on
% cell arrays directly. Also, it removes the issue of drawing diagonals for cells with
% only 2 points of data. This method is based on Jefferey Chiou's raster function.
%
%
% youkitan 2014-11-20
%   ACarey edit (Mar 2015) for colormaps

% Process cfg input
cfg_def.SpikeHeight = 0.4;
cfg_def.axisflag = 'tight';
cfg_def.spkColor = 'k';
cfg_def.lfpColor = 'r';
cfg_def.axislabel = 'on';
cfg_def.openInAxes = 0;
cfg_def.LineWidth = 1;
cfg = ProcessConfig2(cfg_def,cfg_in);

% Check inputs
assert(isfield(S_in,'t'),'Input is not a ts object!')
for i = 1:length(S_in.t)
    assert(size(S_in.t{i},2) == 1,...
        'Poorly constructed ts! Cell %d has incorrect dimensions.',i)
end

%
if cfg.openInAxes % user defined axes
    axes(cfg.openInAxes);
end


spkColor = cfg.spkColor; colorType = 'single color';

% make colour matrix if user specifies a colormap like 'jet' or 'summer'
if ischar(cfg.spkColor) && length(cfg.spkColor) > 1 % if spikeColor is a multicharacter array
    nCells = length(S_in.t);
    spkColor = colormap(eval([cfg.spkColor '(' num2str(nCells) ')'])); % turn the string into a matrix of RGB values
    colorType = 'color matrix';
end

% Plot spikes
[nTrials,nCells] = size(S_in.t);

if  nTrials == 1
    hold on;
    
    for iC = 1:nCells
        nSpikes = length(S_in.t{iC});
        xvals = [ S_in.t{iC}' ; S_in.t{iC}' ; nan(size(S_in.t{iC}')) ];
        yvals = [ iC-cfg.SpikeHeight.*ones(1,nSpikes);...
            iC+cfg.SpikeHeight.*ones(1,nSpikes) ; nan(1,nSpikes) ];
        
        xvals = xvals(:);
        yvals = yvals(:);
        
        switch colorType
            case 'single color'
                plot(xvals,yvals,'Color',spkColor,'LineWidth',cfg.LineWidth)
            case 'color matrix'
                plot(xvals,yvals,'Color',spkColor(iC,:),'LineWidth',cfg.LineWidth)
        end
        
    end
    
    if strcmp(cfg.axislabel,'on')
        % Set axis labels
        ylabel('Cell #','fontsize',10);
        xlabel('Time (sec)','fontsize',10);
    end
    
    % Set ylims
    ylims = [1 nCells];
    
    
else %equivalent to elseif nTrials > 1 (nTrials is always a positive integer)
    cmap = linspecer(nCells);
    hold on;
    for iC = 1:nCells
        for iT = 1:nTrials
            nSpikes = length(S_in.t{iT,iC});
            xvals = [ S_in.t{iT,iC}' ; S_in.t{iT,iC}' ; nan(size(S_in.t{iT,iC}')) ];
            yvals = [ iT-cfg.SpikeHeight.*ones(1,nSpikes);...
                iT+cfg.SpikeHeight.*ones(1,nSpikes) ; nan(1,nSpikes) ];
            
            xvals = xvals(:);
            yvals = yvals(:);
            
            plot(xvals,yvals,'Color',cmap(iC,:),'LineWidth',cfg.LineWidth)
            
        end %iterate trials
    end %iterate cells
    
    if strcmp(cfg.axislabel,'on')
        % Set axis labels
        ylabel('Trial','fontsize',10);
        xlabel('Time (sec)','fontsize',10);
    end
    
    % Set ylims
    ylims = [1 nTrials];
    
end

end