function h = MultiRaster(cfg_in,S)
% function MultiRaster(cfg,S) plots the spiketrain S with optional inputs for lfps and
% iv or ts event objects. 
% User can move viewing window via keyboard input (type "help navigate" in command window).
%
%   INPUTS:
%       cfg_in: input cfg
%       S: input ts NOTE: S.t can be a MxN cell array where M is the number of repeated trials
%       and N is the number of labeled signals
%
%   OUTPUT
%       * handles output depends on plotmode
%        h.S        - handles to spike trains 
%        h.S_iv     - spike train intervals
%        h.LFP      - hande to lfp
%        h.LFP_iv   - handle to LFP intervals
%        h.plotmode - which plot mode was entered
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
%       cfg.axisflag - default 'spandex'
%           Converts axis scaling. TIGHT restricts the plot to the data range. ALL sets
%           the x-axis to cover the whole experiment. 
%           'spandex' tightest possible 
%
%       cfg.LineWidth - default 1
%           Line thickness for each spike.
%
%       cfg.spkColor - default 'k' (black)
%           Color specifier for spikes in a single trial. Can be 1x3 vector of RGB values 
%           [1 0 0], or a string specifying the short name of a color ('r'), or a string 
%           specifying a colormap ('linspecer'). 
%           NOTE: multiple trials will call upon a function to create a good
%           colormap to distinguish cells.
%
%       cfg.ivColor - default 'r' (red)
%           Color specifier for iv. Can be 1x3 vector of RGB values [1 0 0], or a string 
%           specifying the short name of a color ('r').
%
%       cfg.lfpColor - default 'k' (black)
%           Color specifier for lfp signal. Can be 1x3 vector of RGB values [1 0 0], or a string 
%           specifying the short name of a color ('r').
%
%       cfg.lfpHeight - default 5 
%           The height of the lfp from maximum to minium in y-axis units. Vertically 
%           stretches the lfp to improve visibility of oscillations. lfp + iv mode only. 
%
%       cfg.lfpWidth - default 1
%           Line thickness for LFP. Implemented in some modes but not others.
%
%       cfg.lfpMax - default 15
%           Values greater than lfpMax times above mean(abs(lfp.data)) will be plotted 
%           as NaNs. Improves visibilty by cutting off high-amplitude [noise] blips 
%           in raw lfp. lfp + iv mode only.
%
%       cfg.axislabel - default 'on'
%           Automatically creates labels for the axes. Can be turned off.
%      
%       cfg.openNewFig - default 1 
%           If you want to use MultiRaster with subplot, set cfg.openNewFig
%           to 0.Open MultiRaster and read the Help section for instructions.
%           
%
%
% youkitan 2014-11-06 
% edit 2015-01-20
% aacarey edit, 2015-01-20 (lfpHeight and lfpMax)
% aacarey edit Sept 2015, +cfg.openNewFig, removed cfg.openInAxes, added
%      cfg.axisflag option
% youkitan 2016-08-18 multiple lfp fix and ts event fix

%% HELP

% HOW DO I USE MULTIRASTER WITH SUBPLOT? (aacarey)

%     Example:
%     figure('KeyPressFcn',@navigate) % if you want to preserve the x-axis navigation

%     ax1 = subplot(211);
%     MultiRaster(cfg,S1)

%     ax2 = subplot(212);
%     MultiRaster(cfg,S2)

%     linkaxes([ax1 ax2],'x'); % for navigation to work on both axes simultaneously

%% Set cfg parameters and check inputs
cfg_def.SpikeHeight = 0.4;
cfg_def.axisflag = 'spandex';
cfg_def.spkColor = 'k';
cfg_def.LineWidth = 1;
cfg_def.ivColor = 'r';
cfg_def.lfpColor = 'k'; 
cfg_def.lfpHeight = 15;
cfg_def.lfpWidth = 1;
cfg_def.lfpMax = 15;
cfg_def.axislabel = 'on';
cfg_def.windowSize = 1;
cfg_def.openNewFig = 1;
cfg_def.setAxes = 'on';
cfg_def.verbose = 0;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if ~CheckTS(S)
    error('Input spiketrain is not a correctly formed ts.')
end

%% Setup navigate

global evtTimes windowSize time usrfield
windowSize = cfg.windowSize;

% Initialize time vector for navigating
binSize = 0.01;
if ~isfield(cfg,'lfp')
    %create internal tvec that runs from the first spike time from any cell
    %to the last spike time
    tstart = firstSpike(S);
    tend = lastSpike(S);    
    time = tstart:binSize:tend;
else
    for iLFP = 1:length(cfg.lfp)
       if ~CheckTSD(cfg.lfp(iLFP))
           error('Input tsd is not correctly formed.')
       end
    end
    time = cfg.lfp(1).tvec(1):binSize:cfg.lfp(1).tvec(end);
end

% Initialize events
if isfield(cfg,'evt')
    % Initialize usr input text
    if isfield(cfg.evt,'usr')
        usrfield = cfg.evt.usr;
    else
        usrfield = [];
    end
else
    % navigate by seconds
    evtTimes = time(1:100:end); 
end
    
% handle lone figure plotting or subplotting
if cfg.openNewFig
    figure('KeyPressFcn',@navigate)
    hold on;
else 
    hold on;
end


%% Error checking and plot type setup
% Check to see what datatypes we need to plot besides spikes and do error checking. 

ts_only = @(x) isfield(x,'t') && ~isfield(x,'tstart');
iv_only = @(x) isfield(x,'tstart') && ~isfield(x,'t');

if isfield(cfg,'lfp') %lfp
    hS = PlotSpikeRaster2(cfg,S);
    ylims = get(gca,'YLim');
    
    spikelims = get(gca,'Xlim');
    for iLFP = 1:length(cfg.lfp)
        cfg.lfp(iLFP) = lfp_check(cfg,cfg.lfp(iLFP),spikelims);
    end
    
    if ~isfield(cfg,'evt')
        plotMode = 5;
        
    elseif ts_only(cfg.evt) %lfp + ts
        plotMode = 6;
        
        if ~CheckTS(cfg.evt)
            error('cfg.evt is not a correctly formed ts.')
        end
        
    elseif iv_only(cfg.evt) %lfp + iv
        plotMode = 7;
        
        if ~CheckIV(cfg.evt)
            error('cfg.evt is not a correctly formed iv.')
        end
        
    else %lfp + ts + iv
        plotMode = 8;
    end
    
elseif isfield(cfg,'evt')
    if ts_only(cfg.evt) %ts
        plotMode = 2;
        
        if ~CheckTS(cfg.evt)
            error('cfg.evt is not a correctly formed ts.')
        end
        
    elseif iv_only(cfg.evt) %iv
        plotMode = 3;
        
        if ~CheckIV(cfg.evt)
            error('cfg.evt is not a correctly formed iv.')
        end
        
    else %ts + iv
        plotMode = 4;
        
    end
    
else %default (spikes only)
    plotMode = 1;
    
end

%% Choose Plotting mode
switch plotMode        
    case 1 % just spikes
        h.S = PlotSpikeRaster2(cfg,S);
        ylims = get(gca,'YLim');
        
    case 2 % ts data only
        evtTimes = (cfg.evt.t{:});
        h.S = PlotSpikeRaster2(cfg,S);
        PlotTSEvt(cfg,cfg.evt)
        ylims = get(gca,'YLim');
        
    case 3 % iv data only
        evtTimes = (cfg.evt.tstart + cfg.evt.tend)./2;
        S_iv = restrict(S,cfg.evt.tstart,cfg.evt.tend);
        h.S = PlotSpikeRaster2(cfg,S); % plots all spikes
        cfg.spkColor = 'r';
        h.S_iv = PlotSpikeRaster2(cfg,S_iv); % plots event spikes overtop in a different color
        ylims = get(gca,'YLim');
        
    case 4 % ts + iv data NOT WORKING YET
        error('ts + iv event data NOT WORKING YET')
        
    case 5 % lfp data only
        numLFP = length(cfg.lfp);  
        
        if numLFP > 1
            cmap = linspecer(numLFP);
            
            for iLFP = 1:numLFP
                lfp = cfg.lfp(iLFP);
                
                abslfp = abs(lfp.data);
                nans_here = abslfp > cfg.lfpMax*mean(abslfp);
                lfp.data(nans_here) = NaN;
                
                lower_val = -iLFP*cfg.lfpHeight;
                upper_val = lower_val+cfg.lfpHeight;
                
                lfp.data = rescale(lfp.data,lower_val,upper_val);
                h.LFP(iLFP) = plot(lfp.tvec,lfp.data,'Color',cmap(iLFP,:),'LineWidth',cfg.lfpWidth);
            end
        else
            lfp = cfg.lfp;
            
            abslfp = abs(lfp.data);
            nans_here = abslfp > cfg.lfpMax*mean(abslfp);
            lfp.data(nans_here) = NaN;
            
            lower_val = -iLFP*cfg.lfpHeight;
            upper_val = lower_val+cfg.lfpHeight;
            
            lfp.data = rescale(lfp.data,lower_val,upper_val);
            h.LFP = plot(lfp.tvec,lfp.data,'Color',cfg.lfpColor,'LineWidth',cfg.lfpWidth);
        end
        ylims = get(gca,'YLim'); ylims(1) = lower_val;

    case 6 % lfp + ts data
        evtTimes = (cfg.evt.t{:});

        numLFP = length(cfg.lfp);  
        
        if numLFP > 1
            cmap = linspecer(numLFP);
            
            for iLFP = 1:numLFP
                lfp = cfg.lfp(iLFP);
                
                abslfp = abs(lfp.data);
                nans_here = abslfp > cfg.lfpMax*mean(abslfp);
                lfp.data(nans_here) = NaN;
                
                lower_val = -iLFP*cfg.lfpHeight;
                upper_val = lower_val+cfg.lfpHeight;
                
                lfp.data = rescale(lfp.data,lower_val,upper_val);
                h.LFP(iLFP) = plot(lfp.tvec,lfp.data,'Color',cmap(iLFP,:),'LineWidth',cfg.lfpWidth);
            end
        else
            lfp = cfg.lfp;
            
            abslfp = abs(lfp.data);
            nans_here = abslfp > cfg.lfpMax*mean(abslfp);
            lfp.data(nans_here) = NaN;
            
            lower_val = -iLFP*cfg.lfpHeight;
            upper_val = lower_val+cfg.lfpHeight;
            
            lfp.data = rescale(lfp.data,lower_val,upper_val);
            h.LFP = plot(lfp.tvec,lfp.data,'Color',cfg.lfpColor,'LineWidth',cfg.lfpWidth);
        end
        ylims = get(gca,'YLim'); ylims(1) = lower_val;
        PlotTSEvt([],cfg.evt)

    case 7 % lfp + iv data
        evtTimes = (cfg.evt.tstart + cfg.evt.tend)./2;

        cfg_temp.display = 'tsd';
        cfg_temp.bgcol = cfg.lfpColor;
        cfg_temp.fgcol = cfg.ivColor;
        
        numLFP = length(cfg.lfp);  
        
        if numLFP > 1
            cmap = linspecer(numLFP);
            
            for iLFP = 1:numLFP
                lfp = cfg.lfp(iLFP);
                
                abslfp = abs(lfp.data);
                nans_here = abslfp > cfg.lfpMax*mean(abslfp);
                lfp.data(nans_here) = NaN;
                
                lower_val = -iLFP*cfg.lfpHeight;
                upper_val = lower_val+cfg.lfpHeight;
                
                lfp.data = rescale(lfp.data,lower_val,upper_val);
                h.LFP(iLFP) = PlotTSDfromIV(cfg_temp,cfg.evt,lfp);
            end
        else
            lfp = cfg.lfp;
            
            abslfp = abs(lfp.data);
            nans_here = abslfp > cfg.lfpMax*mean(abslfp);
            lfp.data(nans_here) = NaN;
            
            lower_val = -iLFP*cfg.lfpHeight;
            upper_val = lower_val+cfg.lfpHeight;
            
            lfp.data = rescale(lfp.data,lower_val,upper_val);
            h = PlotTSDfromIV(cfg_temp,cfg.evt,lfp);
        end
        ylims = get(gca,'YLim'); ylims(1) = lower_val;
    
    case 8 % lfp + iv + ts data
        error('ts + iv event data NOT WORKING YET')
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
    case 'spandex' % because sometimes tight just isn't tight enough (aacarey sept 2015)
        xlim([time(1) time(end)])
        if isfield(cfg,'lfp');     
            if length(cfg.lfp) > 1
                ylim([-cfg.lfpHeight*length(cfg.lfp) ylims(2)])
            else
                ylim([max(-cfg.lfpHeight) ylims(2)])
            end
        else
            ylim([ylims(1)-1 ylims(2)+1])
        end
end

if strcmp(cfg.setAxes,'off') % (aacarey sept 2015)
    axis off
else
    axis on
end

if cfg.openNewFig; hold off; end


%% output for debugging
plotmodes = {'spikes only','ts events','iv events','ts + iv events','lfp','lfp + ts events'...
    'lfp + iv events','all'};
% out.yvals = yvals;
h.plotMode = plotmodes{plotMode};
% out.plot = h;
if exist('hS','var'); h.S = hS; end
end