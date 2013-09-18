function plot_isi( spikes, show, show_isi )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_isi - plot histogram of ISI distribution for a cluster
%
% Usage:
%       plot_isi( spikes, show, show_isi )
%
% Description:
%   Plots a histogram either representing the inter-spike interval (ISI)
% distribution for the selected spike events or the autocorrelation function
% in Hertz.  Also plotted is a gray area representing the "shadow"
% period during spike detection and a red area representing the user
% defined refractory period.
%
% The width of the histogram bins as well as the maximum time lag displayed
% are set by the following parameters:
%
%    spikes.params.display.isi_bin_size 
%    spikes.params.display.max_isi_to_display
%    spikes.params.display.correlations_bin_size 
%    spikes.params.display.max_autocorr_to_display
%    
% The user can switch between displaying ISIs or autocorrelatoin by right-
% clicking the axes and selecting from a context menu.  The choice can
% also be imposed on all isi plots in the same figure.
%
% On the y-axis is listed the number of refractory period violationss (RPVs)
% along with an estimated contamination percentage and its 95% confidence
% interval under the assumption that contaminating spikes are independent
% events. See poisson_contamination.m for more details.
%
% Inputs:
%   spikes        - a spikes structure
%
% Optional inputs:
%   show          - array describing which events to show in plot
%                 - see get_spike_indices.m, (default = 'all')
%   show_isi      - 1 -> show ISI, 0 -> show autocorrelation
%                 - default is set by spikes.params.display.show_isi
%

    % check arguments
    if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end
    if nargin < 3, show_isi = spikes.params.display.show_isi;  end
    if nargin < 2, show = 'all'; end

    % get the spiketimes
    select = get_spike_indices(spikes, show );      
    spiketimes = sort( spikes.unwrapped_times(select) );
    
    % save data in these axes
    data.spiketimes         = spiketimes;
    data.show_isi           = show_isi;
    data.isi_maxlag         = spikes.params.display.max_isi_to_display;
    data.autocorr_maxlag    = spikes.params.display.max_autocorr_to_display;
    data.shadow             = spikes.params.shadow;
    data.refractory_period  = spikes.params.refractory_period;
    data.corr_bin_size           = spikes.params.display.correlations_bin_size;
    data.isi_bin_size           = spikes.params.display.isi_bin_size;

    % estimate poisson contamination for y-axis
    [expected,lb,ub,RPV] = ss_rpv_contamination( spikes, show  );
    if isempty(lb)
        data.ystr = [num2str(RPV) ' RPVs (' num2str(round(expected*100)) '%)' ];
    else    
        data.ystr = [num2str(RPV) ' RPVs (' num2str(round(lb*100)) '-' num2str(round(expected*100)) '-' num2str(round(ub*100)) '%)' ];
    end
    
    set(gca,'UserData', data,'Tag','isi' )

    % write updating method
    update_isi( [], [], show_isi, gca);
end
    

% callback that updates display on axes
function update_isi( hObject, event, displaymode, ax)

        % get display mode
        data = get( ax,'UserData');
        data.show_isi = displaymode;
        set(ax,'UserData',data)
        set(gcf,'CurrentAxes',ax )
        
        % update display
        make_isi
    
        % set context menu - allow switch between ISI and autocorrelation
        % as well as a global switch for all plot_isi instances on the current figure
        cmenu = uicontextmenu;
        item(1) = uimenu(cmenu, 'Label', 'Show ISI', 'Callback', {@update_isi, 1,ax} );
        item(2) = uimenu(cmenu, 'Label', 'Show autocorr', 'Callback',{@update_isi, 0,ax}  );
        item(3) = uimenu(cmenu, 'Label', 'Use this style on all ISIs in figure', 'Callback',{@impose_all,displaymode,ax},'Separator','on'  );
        set(item(2-displaymode), 'Checked', 'on');    
        set(ax,'UIContextMenu', cmenu )

end

% callback to impose display mode on all plot_isi axes in this figure
function impose_all(hObject,event,displaymode,ax)
        [o,h] = gcbo;
        my_axes = findobj(h,'Tag','isi');
        my_axes = setdiff(my_axes,ax);
        for j = 1:length(my_axes), update_isi([],[],displaymode,my_axes(j)); end
end    
    
% plots the ISI or autocorrelation
function make_isi     

    data = get( gca,'UserData');
    spiketimes = data.spiketimes;
    shadow = data.shadow;
    rp     = data.refractory_period;
    cla reset
    
    % ISI case
    if data.show_isi
        maxlag = data.isi_maxlag;
        bins = round(1000* maxlag/data.isi_bin_size );

        % make plot
        isis = diff(spiketimes);   
        isis = isis(isis <= maxlag); 
        [n,x] =hist(isis*1000,linspace(0,1000*maxlag,bins));
        ymax   = max(n)+1;
        
        % make patches to represent shadow and refractory period

        patch([0 shadow shadow 0 ], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
        patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        hold on,    b2 = bar(x,n,1.0); hold off
        set(b2,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

        % update axes
        set(gca,'YLim',[0 ymax],'XLim',[0 1000*maxlag])
        xlabel('Interspike interval (msec)')
       ylabel({'No. of spikes',data.ystr})
 
    else
        maxlag = data.autocorr_maxlag;
        
        % calculate autocorrelation
        if length(spiketimes) > 1
          [cross,lags] = pxcorr(spiketimes,spiketimes, round(1000/data.corr_bin_size), maxlag);
        else
            cross = 0;  lags = 0;
        end
        cross(find(lags==0)) = 0;
        
        % place patches to represent shadow and refractory period
        ymax = max(cross) + 1;
        patch(shadow*[-1 1 1 -1], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
        patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        patch(-[shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        
        % plot autocorrelation histogram
        hold on, bb = bar(lags*1000,cross,1.0); hold off;  
        set(bb,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

        % set axes
        set(gca, 'XLim', maxlag*1000*[-1 1]);
        set(gca,'YLim',[0 ymax])
        xlabel( 'Time lag (msec)')
        ylabel({'Autocorrelation (Hz)',data.ystr})
    end
    set(gca,'Tag','isi','UserData',data)
    
end