function [ax1,ax2] = plot_stability( spikes, show)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_stability - Plot amplitude and firing rate of a set of spikes over time
%
% Usage:
%      [ax1,ax2] = plot_stability( spikes, clus)
%
% Description:  
%    Creates two axes in the same position to simultaneously plot two 
% aspects of the stability of a given set of spikes.  In the upper half
% and using the right y-axis is a scatter plot of the peak-to-peak amplitude 
% of waveforms taken from the full duration of the recording.  For efficiency
% a maximum number of randomly drawn data points is used as set by
% spikes.params.display.max_scatter.  In the lower half and using the left
% y-axis is a histogram of the firing rate as a function of time.  The 
% unwrapped time is used (see unwrap_time.m).  The size of the bins used
% is set by spikes.params.display.stability_bin_size.
%
% Input:
%   spikes - a spike structure
%
% Optional input:
%   show          - array describing which events to show in plot
%                 - see get_spike_indices.m, (default = 'all')
% 
% Output:
%  ax1 - handle to axes the firing rate over time
%  ax2 - handle to axes that display the amplitude scatter over time
%

   % argument checking
   if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end
   if nargin < 2, show = 'all'; end

   % get data
    select = get_spike_indices(spikes, show );     
    memberwaves = spikes.waveforms(select,:);
    spiketimes  = sort( spikes.unwrapped_times( select ) );
    
    % initialize axes
    cla
    ax1 = gca;
     
    % we are ploting two axes on top of each other
    % the first axes will show the firing rate as a function of time
    
    % bin data
    tlims = [0 sum(spikes.info.detect.dur)];  
    num_bins  = round( diff(tlims) /  spikes.params.display.stability_bin_size);
    edges = linspace(tlims(1),tlims(2),num_bins+1);
    n = histc(spiketimes,edges);  
    n(end) = [];    
    vals = n/mean(diff(edges));
    
    % show histogram
    hold on
    bar(edges(1:end-1) + mean(edges(1:2)),vals,1.0);  shading flat;
    hold off
    set(gca, 'XLim', tlims,'YLim',[0 2*max(get(gca,'YLim'))])
    yticks = get(gca,'YTick'); set(gca,'YTick',yticks( yticks<=max(yticks)/2))
    xlabel('Time (sec)')
    ylabel('Firing rate (Hz)')
    
    % second pass is scatter plot of amplitude versus time
       
    % get amplitudes over time
    amp = range(    memberwaves' );
    
    % only plot a random subset
    if isequal( spikes.params.display.max_scatter, 'all')
        ind = 1:length(amp);
    else
        choice = randperm( length(amp) );
        max_pos = min( length(amp), spikes.params.display.max_scatter );
        ind = choice(1:max_pos );
    end
    
    % prettify axes
    ax2 = axes('Parent',get(ax1,'Parent'),'Unit',get(ax1,'Unit'),'Position',get(ax1,'Position'));
    hold on
    l = scatter(spiketimes(ind),amp(ind));
    hold off
    l1 = line( get(ax2,'XLim'), [ 0 0] );
    set(l1,'Color',[ 0 0 0])
    hold off
    set(l,'Marker','.','MarkerEdgeColor',[.3 .5 .3])
    set(ax2,'Xlim',tlims,'YLim',max( get(ax2,'YLim') ) *[-1 1],'XTickLabel',{},'YAxisLocation','right');
    yticks = get(ax2,'YTick'); set(ax2,'YTick',yticks(yticks>=0))
    ylabel('Amplitude')
    set(ax2,'Color','none')
    
    % linke properties of the 2 axes together
    linkaxes([ax1 ax2],'x');
    linkprop([ax1 ax2],'Position');
    linkprop([ax1 ax2],'Unit');
    
