function out = MultiPlotSpikeRaster(cfg_in,S)
% function MultiPlotSpikeRaster(cfg,S)
%
%   INPUTS:
%       cfg_in: input cfg
%       S: input ts NOTE: S.t can be a MxN cell array where M is the number of repeated trials
%       and N is the number of labeled signals
%
%   CFG OPTIONS:
%       cfg.SpikeHeight - default 0.4
%           Size of marker for each spike
%
%       cfg.lfp - default NONE
%           LFP data in ts format. Multiple sources of lfp data can be input as a single
%           cfg.lfp object by using a 1xM struct array where M is the number of lfp signals.
%           (i.e., cfg.lfp(1) = lfpA; cfg.lfp(2) = lfpB)
%
%       cfg.evt - default NONE
%           Event data in ts or iv format. Multiple sources of event data in a single 
%           cfg.evt object by using a 1xM cell array where M is the number of event files.
%           (i.e., cfg.evt = {evt1,evt2}, where evt1=ts evt2=iv)
%
%       cfg.axisflag - default 'tight'
%           Converts axis scaling. TIGHT restricts the plot to the data range. ALL sets
%           the x-axis to cover the whole experiment.
%
%       cfg.spkColor - default 'k' (black)
%           Color specifier for spikes in a single trial. Can be any MATLAB string color
%           specifier. NOTE: multiple trials will call upon a function to create a good
%           colormap to distinguish cells.
%
%       cfg.lfpColor - default 'r' (red)
%           Color specifier for lfp signal. Can be any MATLAB string color specifier.
%
%       cfg.axislabel - default 'on'
%           Automatically creates labels for the axes. Can be turned off.
%
%
% youkitan 2014-11-06 

%% Set cfg parameters and check inputs
cfg_def.SpikeHeight = 0.4;
cfg_def.axisflag = 'tight';
cfg_def.spkColor = 'k';
cfg_def.lfpColor = 'r';
cfg_def.axislabel = 'on';
cfg = ProcessConfig2(cfg_def,cfg_in);

%% Extra plotting
% Check to see what datatypes we need to plot besides spikes and do error checking. 

ts_only = @(x) isfield(x,'t') && ~isfield(x,'tstart');
iv_only = @(x) isfield(x,'tstart') && ~isfield(x,'t');

if isfield(cfg,'lfp') %lfp
    ylims = PlotSpikeRaster2(cfg,S);

    spikelims = get(gca,'Xlim');
    for iLFP = 1:length(cfg.lfp)
      cfg.lfp(iLFP) = lfp_check(cfg,cfg.lfp(iLFP),spikelims);
    end

    if ~isfield(cfg,'evt')
        plotMode = 5;

    elseif ts_only(cfg.evt) %lfp + ts
        plotMode = 6;

    elseif iv_only(cfg.evt) %lfp + iv
        plotMode = 7;

    else %lfp + ts + iv
        plotMode = 8;
    end

elseif isfield(cfg,'evt')
    if ts_only(cfg.evt) %ts
        plotMode = 2;

    elseif iv_only(cfg.evt) %iv
        plotMode = 3;

    else %ts + iv
        plotMode = 4;

    end

else %default (spikes only)
    plotMode = 1;

end


%% Plotting extra data
switch plotMode        
    case 1 % just spikes
        ylims = PlotSpikeRaster2(cfg,S);
        
    case 2 % ts data only
        ylims = PlotSpikeRaster2(cfg,S);
        PlotTSEvt(cfg,cfg.evt)
        
    case 3 % iv data only
        S_iv = restrict(S,cfg.evt.tstart,cfg.evt.tend);
        ylims = PlotSpikeRaster2(cfg,S);
        cfg.spkColor = 'r';
        ylims = PlotSpikeRaster2(cfg,S_iv);
                    
    case 4 % ts + iv data NOT WORKING YET
        for i = 1:length(cfg.evt)
        end
        
        S_iv = restrict(S,cfg.evt.tstart,cfg.evt.tend);
        PlotSpikeRaster2(cfg,S);
        cfg.spkColor = 'r';
        [xvals,yvals,ylims] = PlotSpikeRaster2(cfg,S_iv);
        PlotTSEvt(cfg.evt,ylims)
        
    case 5 % lfp data only
        numLFP = length(cfg.lfp);        
        cmap = linspecer(numLFP);

        for iLFP = 1:numLFP
            lfp = cfg.lfp(iLFP);
            lower_val = -iLFP*5;
            upper_val = lower_val+5;

            lfp.data = rescale(lfp.data,lower_val,upper_val);
            plot(lfp.tvec,lfp.data,'Color',cmap(iLFP,:))
        end
        ylims(1) = lower_val;

    case 6 % lfp + ts data
        cfg.lfp.data = rescale(cfg.lfp.data,-10,0);
        plot(cfg.lfp.tvec,cfg.lfp.data,'Color',cfg.lfpColor);
        ylims(1) = -10;
        PlotTSEvt([],cfg.evt)

    case 7 % lfp + iv data
        cfg.lfp.data = rescale(cfg.lfp.data,-10,0);
        PlotTSDfromIV([],cfg.evt,cfg.lfp);
        ylims(1) = -10;
    
    case 8 % lfp + iv + ts data
        [xvals,yvals,ylims] = PlotSpikeRaster2(cfg,S);

end 

%% Helper Functions

    function lfp = lfp_check(cfg_in,lfp_in,limits)
    % Checks lfp size and restricts it to the spiking data time interval
        if (min(lfp_in.tvec) < limits(1) || max(lfp_in.tvec) > limits(2))...
                && strcmp(cfg_in.axisflag,'tight')
            fprintf('minlfp = %d maxlfp = %d minspikes = %d maxspikes = %d\n',...
                min(lfp_in.tvec),max(lfp_in.tvec),limits(1),limits(2))
            fprintf('Range of lfp data exceeds spiketrain data!... Restricting lfp data...');
            lfp = restrict(lfp_in,limits(1),limits(2));   
        else 
            lfp = lfp_in;
        end
    end


%% Adjust figure

% Set Axis limits
xlims = get(gca,'XLim');

switch cfg.axisflag
    case 'tight'
        xlim([xlims(1) xlims(2)]);
        ylim([ylims(1)-1 ylims(2)+1]);
    case 'all'
        xlim([S.cfg.ExpKeys.TimeOnTrack S.cfg.ExpKeys.TimeOffTrack]);
        ylim([ylims(1)-1 ylims(2)+1])
end

hold off;



%% output for debugging
plotmodes = {'spikes only','ts events','iv events','ts + iv events','lfp','lfp + ts events'...
    'lfp + iv events','all'};
% out.yvals = yvals;
out.plotMode = plotmodes{plotMode};
% out.plot = h;
end