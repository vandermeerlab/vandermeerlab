function PlotSpikeRaster(cfg_in,S)
% function PlotSpikeRaster(cfg,S)
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
%           lfp data as a tsd
%
%       cfg.axisflag - default 'tight'
%           Converts axis scaling. TIGHT restricts the plot to the data range. ALL sets
%           the x-axis to cover the whole experiment.
%
%       cfg.Color - default 'k' (black)
%           Color specifier for spikes in a single trial. Can be any MATLAB string color
%           specifier. NOTE: multiple trials will call upon a function to create a good
%           colormap to distinguish cells.
%
%
% MvdM 2014-07-20
% youkitan 2014-11-06 edit

%% Set cfg parameters
cfg_def.SpikeHeight = 0.4;
cfg_def.axisflag = 'tight';
cfg_def.Color = 'k';
cfg = ProcessConfig2(cfg_def,cfg_in);


%% Plot spikes

%Check number of runs (trials) and cells
[nTrials,nCells] = size(S.t);

if  nTrials == 1
    hold on;
    for iC = 1:nCells
        if isempty(S.t{iC})
            continue;
        end

        h{iC} = plot([S.t{iC} S.t{iC}],[iC-cfg.SpikeHeight iC+cfg.SpikeHeight],cfg.Color);
    end
    
    % Set axis labels
    ylabel('Cell #','fontsize',12);
    xlabel('Time (sec)','fontsize',12);
    
    % Set ylims
    ylims = [1 nCells];

    
else %equivalent to elseif nTrials > 1 (nTrials is always a positive integer)
    cmap = linspecer(nCells);
    hold on;
    for iC = 1:nCells
        for iT = 1:nTrials
            plot([S.t{iT,iC} S.t{iT,iC}],[iT-cfg.SpikeHeight iT+cfg.SpikeHeight],'Color',cmap(iC,:))
        end %iterate trials
    end %iterate cells
    
    % Set axis labels
    ylabel('Trial','fontsize',12);
    xlabel('Time (sec)','fontsize',12);
    
    % Set ylims
    ylims = [1 nTrials];

end


%% Plot LFP

if isfield(cfg,'lfp')
    
    spikelims = get(gca,'Xlim');
    
    % Make sure lfp data has same time interval as spiking data
    if (min(cfg.lfp.tvec) < spikelims(1) || max(cfg.lfp.tvec) > spikelims(2))...
            && strcmp(cfg.axisflag,'tight')
        fprintf('minlfp = %d maxlfp = %d minspikes = %d maxspikes = %d\n',...
            min(cfg.lfp.tvec),max(cfg.lfp.tvec),spikelims(1),spikelims(2))
        fprintf('Range of lfp data exceeds spiketrain data!... Restricting lfp data...');
        cfg.lfp = restrict(cfg.lfp,spikelims(1),spikelims(2));
    end

    cfg.lfp.data = rescale(cfg.lfp.data,-5,0);
    plot(cfg.lfp.tvec,cfg.lfp.data,'Color',[0.5 0.5 0.5]);
end


%% Adjust figure

% Set Axis limits
xlims = get(gca,'Xlim');

switch cfg.axisflag
    case 'tight'
        xlim([xlims(1) xlims(2)]);
        ylim([ylims(1)-1 ylims(2)+1]);
    case 'all'
        xlim([S.cfg.ExpKeys.TimeOnTrack S.cfg.ExpKeys.TimeOffTrack]);
        ylim([ylims(1)-1 ylims(2)+1])
end

hold off;

end