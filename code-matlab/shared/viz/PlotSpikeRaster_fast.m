function h = PlotSpikeRaster2(cfg_in,S_in)
%% Plot spikes
% Converting the data to (x,y) coordinates is significantly faster to plot than calling on
% cell arrays directly. Also, it removes the issue of drawing diagonals for cells with
% only 2 points of data. This method is based on Jefferey Chiou's raster function.
%
%
% MvdM adapted from youkitan & aacarey PlotSpikeRaster()

% Process cfg input
cfg_def.axisflag = 'tight';
cfg_def.spkColor = 'k';
cfg_def.lfpColor = 'r';
cfg_def.axislabel = 'on';
cfg_def.openInAxes = 0;
cfg_def.LineWidth = 1;
cfg_def.verbose = 0;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

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

% handle color matrix (aacarey)
if length(cfg.spkColor) == 1 || isnumeric(cfg.spkColor)
    cfg.spkColor = repmat(cfg.spkColor,length(S_in.t),1);
    
elseif length(cfg.spkColor) > 1 && ~isnumeric(cfg.spkColor)
    cfg.spkColor = colormap(eval([cfg.spkColor '(' num2str(length(S_in.t)) ')'])); % turn the string into a matrix of RGB values
else
    error('Unable to assign cfg.spkColor')
end

% Plot spikes
[nTrials,nCells] = size(S_in.t);


hold on;
%     h = nan(size(S_in.t));
h = cell(size(S_in.t));
for iC = 1:nCells
    nSpikes = length(S_in.t{iC});
    xvals = S_in.t{iC}';
    yvals = iC.*ones(1,nSpikes);
    
    xvals = xvals(:);
    yvals = yvals(:);
    
    if isempty(xvals) && isempty(yvals)
        %             h(iC) = nan;
        h{iC} = nan;
    else
        %             h(iC) = plot(xvals,yvals,'Color',cfg.spkColor(iC,:),'LineWidth',cfg.LineWidth);
        plot(xvals,yvals,'.','Color',cfg.spkColor(iC,:));
        h{iC} = get(gca);
    end
end

if strcmp(cfg.axislabel,'on')
    % Set axis labels
    ylabel('Cell #','fontsize',10);
    xlabel('Time (sec)','fontsize',10);
end

% Set ylims
ylims = [0.5 nCells+0.5];



set(gca,'YLim',ylims)

end