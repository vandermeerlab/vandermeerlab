function [cross,lags,collision_info] = plot_xcorr( spikes, show1, show2 )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_xcorr - plot cross-correlation between 2 spike clusters
%
% Usage:
%       [cross,lags,collision_info] = plot_xcorr( spikes, show1, show2 )
%
% Description:
%   Plots a histogram either representing the cross-correlation between
% the two spike trains of the specified clusters or representing the 
% autocorrelation function if the 2 spike trains are merged, either are in 
% Hz.  Also plotted is a gray area representing the "shadow"
% period during spike detection and a red area representing the user
% defined refractory period.
%
% The width of the histogram bins as well as the maximum time lag displayed
% are set by the following parameters:
%
%    spikes.params.display.correlations_bin_size 
%    spikes.params.display.max_autocorr_to_display
%    
% The user can switch between displaying ISIs or autocorrelatoin by right-
% clicking the axes and selecting from a context menu.  The choice can
% also be imposed on all isi plots in the same figure.  The default display
% mode is set by:
%
%     spikes.params.display.default_xcorr_mode
%
% On the y-axis is listed the number spikes in both clusters.  Along the
% title are 3 staticis about refractory period violations (RPVs0.
%
%       TOT - RPVs when both spike trains are merged 
%       NEW - No. of RPVs created when the spike trains are merged (above what is already present)
%       EXP - No. of RPVs expected if 2 spike trains are independent
%
% Inputs:
%   spikes        - a spikes structure
%   show1         - spikes selected for group 1 (see get_spike_indices.m) 
%   show2         - spikes selected for group 2
%
% Output:
%   cross             - cross-correlation between show1 and show2 (Hz)
%   lags              - time lags associated with cross (s)
%   collision_info    - structure with 3 fields (tot,new,exp) as defined above
%

   % check parameters
   if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end
   select1 = get_spike_indices(spikes, show1 );
   select2 = get_spike_indices(spikes, show2 );

    % set internal data
    data.st1 = spikes.unwrapped_times(select1);
    data.st2 = spikes.unwrapped_times(select2);
    data.show_xcorr = spikes.params.display.default_xcorr_mode;
    data.shadow             = spikes.params.shadow;
    data.refractory_period  = spikes.params.refractory_period;
    data.corr_bin_size           = spikes.params.display.correlations_bin_size;
    data.maxlag = spikes.params.display.max_autocorr_to_display;  % (msec);
    
    % collect RPV statistics
    RPV1 = sum(diff(data.st1) < data.refractory_period / 1000 );
    RPV2 = sum(diff(data.st2) < data.refractory_period / 1000 );
    st12 = sort( [data.st1 data.st2] );
    data.tot_rpv = sum(diff(st12) < data.refractory_period / 1000 );
    data.new_rpv = data.tot_rpv - (RPV1+RPV2);
    data.exp_rpv = round( 2*((data.refractory_period - data.shadow)/1000)*length(data.st1)*length(data.st2)/sum(spikes.info.detect.dur) );
    
    collision_info.tot = data.tot_rpv;
    collision_info.new = data.new_rpv;
    collision_info.exp = data.exp_rpv;
    
    % prepare axes
    cla   
    set(gca,'UserData', data,'Tag','xcorr' )

    % call update callback
    update_xcorr( [], [], data.show_xcorr, gca);
end 
    
% callback for updating axes when user selects new display mode
function update_xcorr( hObject, event, displaymode, ax)

    % set up
    data = get( ax,'UserData');
    data.show_xcorr = displaymode;
    set(ax,'UserData',data)
    set(gcf,'CurrentAxes',ax )

    % update display
    make_xcorr
    
    % set context menu
    cmenu = uicontextmenu;
    item(1) = uimenu(cmenu, 'Label', 'Show x-corr of spike trains', 'Callback', {@update_xcorr, 1,ax} );
    item(2) = uimenu(cmenu, 'Label', 'Show autocorr of merged spike trains', 'Callback',{@update_xcorr, 0,ax}  );
    item(3) = uimenu(cmenu, 'Label', 'Use this style on all ISIs in figure', 'Callback',{@impose_all,displaymode,ax},'Separator','on'  );
    
    set(item(2-displaymode), 'Checked', 'on');    
    set(ax,'UIContextMenu', cmenu )
end

% callback to impose display mode on all xcorr plots in figure
function impose_all(hObject,event,displaymode,ax)
        [o,h] = gcbo;
        my_axes = findobj(h,'Tag','xcorr');
        my_axes = setdiff(my_axes,ax);
        for j = 1:length(my_axes), update_xcorr([],[],displaymode,my_axes(j)); end
end
    
% updates the display
function make_xcorr     

    % get some numbers
    data = get( gca,'UserData');
    shadow = data.shadow;
    rp     = data.refractory_period;
    maxlag = data.maxlag;

    cla
    
    %  choose cross-correlation or merged autocorrelation
    if data.show_xcorr
        st1 = data.st1;
        st2 = data.st2;
        ystr{1} = 'Cross-corr (Hz)';
    else
        st1 = sort( [data.st1 data.st2] );
        st2 = st1;
        ystr{1} = 'Merged Autocorr (Hz)';
    end

    ystr{2} = ['N_1 = ' num2str(length(data.st1)) ', N_2 = ' num2str(length(data.st2))];

    if length(st1) > 1 & length(st2) > 1
      [cross,lags] = pxcorr(st1,st2, round(1000/data.corr_bin_size), maxlag);
    else
        cross = 0;  lags = 0;
    end

    cross(find(lags==0)) = 0;

    % place  shadow and refractory period areas
    ymax = max(cross) + 1;
    patch(shadow*[-1 1 1 -1], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
    patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
    patch(-[shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
    
    % place xcorr histogram
    hold on, bb = bar(lags*1000,cross,1.0); hold off;  
    set(bb,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

    % update axes
    set(gca, 'XLim', maxlag*1000*[-1 1]);
    set(gca,'YLim',[0 ymax])
    xlabel( 'Time lag (msec)')
    ylabel(ystr)
    title( [ 'RPVs: ' num2str(data.tot_rpv) ' TOT, ' num2str(data.new_rpv) ' NEW, '  num2str( data.exp_rpv ) ' EXP']);  
    set(gca,'Tag','xcorr','UserData',data)

end