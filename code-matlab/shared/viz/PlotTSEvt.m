    function [] = PlotTSEvt(cfg_in,evt_in)
    %% PLOTTSEVT Plot TS Events
    %   Draws vertical, color coded lines for timestamped events. Use this to plot events on
    %   top of a existing spike raster.
    %
    %   INPUT:
    %       cfg_in: input cfg
    %       evt_in: input ts
    %   
    %   EXAMPLE:
    %       S = LoadSpikes(cfg)
    %       evt = LoadEvents(cfg)
    %       plotspikes(cfg,S)
    %       PlotTSEvt(cfg,evt)
    %
    %
    % youkitan 2014-11-20
    
    %% Set cfg parameters and check inputs
    cfg_def = [];
    cfg_def.legend = 'off';
    cfg_def.openInAxes = 0; % may not need this, but is more transparent?
    cfg = ProcessConfig2(cfg_def,cfg_in);
    
    %% Plot
    events = evt_in.t;
    labels = evt_in.label;
    numEvt = length(events);
    cmap = linspecer(numEvt);
    h = [];
    legendnames = {};

    if ~isempty(findall(0,'Type','Figure')) %checks to see if there is a figure open
        ylims = get(gca,'YLim');
    else
        ylims = [0 10];
    end

    hold on;
    for iEvt = 1:length(events)                     
        numEvt = length(events{iEvt});

        xvals = [ events{iEvt} ; events{iEvt} ; nan(size(events{iEvt})) ];
        yvals = [ ylims(1).*ones(1,numEvt); ylims(2).*ones(1,numEvt) ; nan(1,numEvt) ];

        xvals = xvals(:);
        yvals = yvals(:);

        h(iEvt) = plot(xvals,yvals,':','LineWidth',1,'Color',cmap(iEvt,:));
        legendnames{iEvt} = [labels{iEvt}];

    end %iterate events
    
    if strcmp(cfg.legend,'on')
        legend(h,legendnames,'Interpreter','none','Location','BestOutside')
    end
    
    end