function plot_features(spikes, show, colormode, show_outliers, alt_assigns )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_features - interactive 2D scatter plot of spike features
%
% Usage:
%      plot_features(spikes, show, colormode, show_outliers, alt_assigns )
%
% Description:  
%     This function plots a 2D scatter plot of various features of spike events. 
% There are several features available and some of them require additional
% arguments.
%
%  Feature list:
%      Signal        - value of spike waveform at a particular sample
%                    - The sample number is a required argument.  If the waveform has multiple 
%                      electrodes, then the sample number is for the concatenated waveform.
%      PC            - value of principal componenet
%                    - The principal component number is a required argument.
%      Cluster       - cluster ID
%      Minicluster   - minicluster ID 
%      Event Time    - "unwrapped" time of event
%      ISI Preceding - inter-spike interval (ms) between current event and previous ones of the same cluster
%                    - in order to focus on short latencies, time is truncated at spikes.params.display.max_isi_to_display 
%      Total Energy  - norm of the waveform (i.e., square root of sum of squares)
%      Amplitude     - max value minus min value of waveform
%      Width         - time difference in samples between maximum and minimum of waveform
%
%  In general, if cluster IDs are not available, the feature will use minicluster IDs.  
%  If minicluster IDs are not available, all events will use same ID.  The default 
%  feature for the x and y axes and the argument to these features are set by:
%
%      spikes.params.display.xchoice         
%      spikes.params.display.xparam  
%      spikes.params.display.ychoice  
%      spikes.params.display.yparam
%
%  The currently displayed feature can be changed via a pop-up menu that you can access
%  by clicking on the axis label.  Note that if the default feature is unavailable, such as
%  with the choice of "Cluster" or "Minicluster", then the default is to use "PC" 1 or 2.
%
%  To change the color mode, show/hide the legend, or show/hide outliers, right-click
%  within the plot to bring up a context-menu.  You will also have the option to 
%  replace the scatter plot with a head map.  The specific color map used is taken from
%  spikes.params.display.cmap.
%
%  To bring a cluster of data points to the front, left click one of its data points.
%  To send a cluster of data points  to the back, right click one of its data points.
%
% Input:
%   spikes - a spike structure
%
% Optional input:
%   show          - array describing which events to show in plot
%                 - see get_spike_indices.m, (default = 'all')
%   colormode     - mode for coloring individual data points 
%                 - 1 => all black, 2 => use color of mincluster, 3 => use color of cluster  
%                 - the default is to use the highest flag that is valid for the given spikes object
%   show_outliers - whether to display events marked as outliers (default is set by spikes.params.display.show_outliers) 
%   alt_assigns   - The cluster ID assignments for each unit (default is to use spikes.assigns)
%                 - If minicluster assignments are present in the spikes object, then alt_assigns
%                   needs to be consistent with this. Particularly, there should not be more IDs
%                   in alt_assigns than in spikes.info.kmeans.assigns.
%                 - If minicluster assignments are absent in the spikes object, then alt_assigns 
%                    will also be used as the minicluster IDs.  
%

    % argument checking
    if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end
    has_miniclusters = nargin ==5 | isfield(spikes.info,'kmeans');
    has_clusters     = nargin ==5 | isfield(spikes,     'assigns');        
    valid_modes = [1 has_miniclusters has_clusters];

    if nargin == 5,  % if given alternate assignments, use them
        spikes.assigns = alt_assigns;
        if ~isfield( spikes.info,'kmeans' ), 
            spikes.info.kmeans.assigns = alt_assigns;
            spikes.info.kmeans.colors = emphasize_cluster_colors(spikes);
        end
    end
    if nargin < 4 , show_outliers  = spikes.params.display.show_outliers;  end      
    if nargin < 3 , colormode      = find(valid_modes,1,'last' );     end
    if nargin < 2,  show           = 'all';                                end
    
    % make this plot's internal data structure    
    show = get_spike_indices(spikes, show );
    data.colormode = colormode;
    data.valid_modes = valid_modes;
    data.show = show;
    data.show_outliers = show_outliers;
    data.show_density = 0;
    data.cmap = spikes.params.display.cmap;
    
    % calculate all the ISI's now because spikes may be separated from their cluster mates later
    spikes = add_isis_field( spikes );
        
    % only store spikes that were specified to be shown
    if ~isequal(show,'all')
      spikes.spiketimes = spikes.spiketimes(show); 
      spikes.unwrapped_times = spikes.unwrapped_times(show); 
      spikes.trials = spikes.trials(show);
      spikes.isis   = spikes.isis(show);
      spikes.waveforms = spikes.waveforms(show,:,:);
      if has_miniclusters, spikes.info.kmeans.assigns = spikes.info.kmeans.assigns(show); end
      if has_clusters, spikes.assigns = spikes.assigns(show); end
    end
    data.spikes = spikes;

    % set feature controls
    xc = spikes.params.display.xchoice; xp = spikes.params.display.xparam;
    yc = spikes.params.display.ychoice; yp = spikes.params.display.yparam;
    if ( ~isfield(spikes.info, 'kmeans' ) & ( isequal(xc, 'Minicluster') ) ) | ( ~isfield(spikes, 'assigns' ) & ( isequal(xc, 'Cluster') ) )
        xc = 'PC'; xp = 1;
    end
    if ( ~isfield(spikes.info, 'kmeans' ) & ( isequal(yc, 'Minicluster') ) ) | ( ~isfield(spikes, 'assigns' ) & ( isequal(yc, 'Cluster') ) )
        yc = 'PC'; yp = 2;
    end
    
    data.xchoice = xc;
    data.xparam = xp;
    data.xcontrol = [];
    data.ychoice = yc;
    data.yparam = yp;
    data.ycontrol = [];
    data.show_legend = 0 ;

    %  initialize axes, figure, and plot handles
    hax = gca;
    cla
    set(hax,'UserData',data);
    hfig = get(hax,'Parent');
    hold on
    if has_miniclusters
        
        miniclusters = sort( unique(spikes.info.kmeans.assigns) );
       
        % make handles for mini-clusters
        for j = 1:length(miniclusters)    
             num_spikes = sum( spikes.info.kmeans.assigns == miniclusters(j) ); 
             h  = plot(zeros([1 num_spikes]), zeros([1 num_spikes]), '.' );
             set(h,'UserData',miniclusters(j), 'ButtonDownFcn', {@raise_me});
        end
        
        % make handle for outliers
        if isfield( spikes.info,'outliers')
             num_spikes = length(spikes.info.outliers.spiketimes); 
             h  = plot(zeros([1 num_spikes]), zeros([1 num_spikes]), '.' );
             set(h,'UserData',0, 'ButtonDownFcn', {@raise_me});
        end
    else
          % only make 1 handle if there are no minicluster
          num_spikes = length(spikes.spiketimes);
          h  = plot(zeros([1 num_spikes]), zeros([1 num_spikes]), '.' );
    end
    
    % for 2-D densities, we maintain 2 separate images, with and without ouliers
    set( image([0 0;0 0]),'Visible','off','Tag','image_with_outliers');
    set( image([0 0;0 0]),'Visible','off','Tag','image_without_outliers');
    
    hold off      
    xlabel('temp'); ylabel('temp');
    
    % now populate the plot
    redraw_features(hax);
    recolor_features(hax);
    
    % set up the call backs
    set(get(hax,'Xlabel'), 'ButtonDownFcn', {@make_control, 'x'});
    set(get(hax,'Ylabel'), 'ButtonDownFcn', {@make_control, 'y'});
    set(hax,'ButtonDownFcn', @delete_controls, 'DeleteFcn', @delete_controls);
    set(hax,'Tag','feature_plot');
    set(hax,'UIContextMenu',  get_menu( data.colormode, data.valid_modes, show_outliers, data.show_legend, data.show_density ) )

end

% callback for when user changes the color mode
function update_colormode( varargin )
    
    % grab a few things
    [o,h] = gcbo;
    hax = findobj(h,'Tag','feature_plot');
    data = get(hax,'UserData');
    data.colormode = varargin{3};
    set(hax,'UserData', data);
   
    % update menu checklist
    cm = get(hax,'UIContextMenu');
    checked = {'off', 'off', 'off'};
    checked{data.colormode} = 'on';
    for j = 1:3, set( findobj(cm,'Tag',['coloroption' num2str(j)]), 'Checked', checked{j}); end 
   
    recolor_features(hax);
end

% helper function to set up context menu
function cmenu = get_menu( colormode, valid_modes, show_outliers, show_legend, show_density )
    checked = { 'off', 'off', 'off' }; checked{colormode} = 'on';
    enabling = { 'on', 'on', 'off' };  if valid_modes(3), enabling{3} = 'on'; end

    has_assigns = any(valid_modes(2:3));
    if show_outliers, out_check = 'on'; else, out_check = 'off'; end
    if show_legend, show_leg ='on'; else, show_leg = 'off';end
    if show_density, density_check ='on'; else, density_check = 'off';end
    
    cmenu = uicontextmenu;
    if has_assigns
      uimenu(cmenu, 'Label', 'Do not color', 'Callback', {@update_colormode, 1} ,'Enable', enabling{1} ,'Checked',checked{1},'Tag','coloroption1');
      uimenu(cmenu, 'Label', 'Color by miniclusters', 'Callback',{@update_colormode, 2},'Enable',enabling{2},'Checked',checked{2},'Tag','coloroption2');
      uimenu(cmenu, 'Label', 'Color by clusters', 'Callback', {@update_colormode, 3},'Enable',enabling{3},'Checked',checked{3},'Tag','coloroption3');
      uimenu(cmenu, 'Label', 'Show legend', 'Callback', @toggle_legend,'Separator','on','Checked',show_leg,'Tag','legendoption');
      uimenu(cmenu, 'Label', 'Show outliers', 'Callback',@toggle_show_outliers,'Checked',out_check,'Tag','showoutlieroption'  );
    end
    uimenu(cmenu, 'Label', 'Show density', 'Callback',@toggle_show_density,'Checked',density_check,'Tag','showdensityoption'  );
    
end

% callback for when user toggles show/hide outliers
function toggle_show_outliers(varargin)

    % change flag
    hax = gca;
    figdata = get( hax,'UserData');
    figdata.show_outliers = ~figdata.show_outliers;
    set( hax,'UserData',figdata);

    % update checkmark
    if figdata.show_outliers, out_check = 'on'; else, out_check = 'off'; end
   
    set( findobj(get(hax,'UIContextMenu'),'Tag','showoutlieroption'), 'Checked', out_check);
    
    % change which picture is visible
    set( findobj(hax,'Tag','image_with_outliers'), 'Visible', bool_to_word( figdata.show_outliers & figdata.show_density ) );
    set( findobj(hax,'Tag','image_without_outliers'), 'Visible', bool_to_word( ~figdata.show_outliers & figdata.show_density ) );
        
    % redraw
    recolor_features(hax)
    axis auto
end

% callback for when user toggles show/hide density
function toggle_show_density( varargin )
 % change flag
    hax = gca;
    figdata = get( hax,'UserData');
    figdata.show_density = ~figdata.show_density;
    set( hax,'UserData',figdata);

    % update checkmark
    out_check = bool_to_word( figdata.show_density );
   
    set( findobj(get(hax,'UIContextMenu'),'Tag','showdensityoption'), 'Checked', out_check);
    
    % redraw image, with or without oultiers
     hndls = findobj(hax,'Type','line');
     if figdata.show_outliers
         item = findobj(hax,'Tag','image_with_outliers');
     else
         item = findobj(hax,'Tag','image_without_outliers');
     end
     
     set( item, 'Visible', bool_to_word( figdata.show_density) );
     if figdata.show_density
         map = figdata.cmap;
         colormap(map); 
         set(gca,'Color',map(1,:))
     else
        recolor_features(hax)
        set(gca,'Color',[ 1 1 1])
     end
     
     % hide legend
      legend hide
      figdata.show_legend = 0; set(hax,'UserData',figdata);
      set( findobj(get(hax,'UIContextMenu'),'Tag','legendoption'), 'Checked', 'off');
      axis auto
end

% utility to turn 1 into 'on' and 0 into 'off'
function word = bool_to_word( b )
    if b, word = 'on'; else, word = 'off'; end
end

% callback for when user toggles show/hide legend
function toggle_legend(varargin)

    % change flag
    hax = gca;
    axdata = get( hax,'UserData');
    axdata.show_legend = ~axdata.show_legend;
    set( hax,'UserData',axdata);

    % update checkmark
    if axdata.show_legend 
        legend('show');  show_leg='on';
    else
        legend('hide'); show_leg='off';
    end
   
    set( findobj(get(hax,'UIContextMenu'),'Tag','legendoption'), 'Checked', show_leg);
    
end
    
% draw the features, one handle at a time
function redraw_features(hax)

    axes(hax)
    data = get(hax,'UserData');    
    hndls = findobj(hax,'Type','line');
    outlier_loc = [];
    
    % update the data plots
    for j = 1:length(hndls)    
        clus = get( hndls(j),'UserData');
            if isempty(clus)
                indices = 1:length(data.spikes.spiketimes);
            elseif clus == 0
                indices = 0;
                outlier_loc = [outlier_loc j];
            else
              indices = find( data.spikes.info.kmeans.assigns == clus );
            end
            x{j} = get_feature( data.spikes,  indices, data.xchoice, data.xparam );
            y{j} = get_feature( data.spikes,  indices, data.ychoice, data.yparam );        
            set( hndls(j),'XData',x{j},'YData',y{j});           
    end
    
        % make the density image, with outliers
        x2 = concatenate_from_cell(x);
        y2 = concatenate_from_cell(y);
        D  = 2*round( min( max( sqrt( length(x2) ) / 2, 10 ), 100 ) );
        [counts,x_inds,y_inds] = histxy(x2, y2, D,2);
        set( findobj(hax,'Tag','image_with_outliers'), 'CData',counts','XData',x_inds,'YData',y_inds)

        % make the density image, without outliers
        keepers = setdiff(1:length(hndls),outlier_loc);
        x2 = concatenate_from_cell(x(keepers));
        y2 = concatenate_from_cell(y(keepers));
        D  = 2*round( min( max( sqrt( length(x2) ) / 2, 10 ), 100 ) );
        [counts,x_inds,y_inds] = histxy(x2, y2, D,2);
        set( findobj(hax,'Tag','image_without_outliers'), 'CData',counts','XData',x_inds,'YData',y_inds)

        set( get(hax,'XLabel'), 'String', [data.xchoice ' ' num2str(data.xparam)]);
        set( get(hax,'YLabel'), 'String', [data.ychoice ' ' num2str(data.yparam)]);
    
end

% helper function to change my cell array into a matrix
function y = concatenate_from_cell( x )
    y = [];
    for j = 1:length(x)
        y = [ y x{j}(:)'];
    end
end

% recolor features, 1 handle at a time
function recolor_features(hax)
    
    axes(hax)
    data = get(hax,'UserData');    
    hndls = findobj(hax,'Type','line');
    
    % for each handle    
    for j = 1:length(hndls)    
        clus = get( hndls(j),'UserData');
        
        % we are in the outlier cluster, paint it black
        if clus == 0
            if data.show_outliers, visible = 'on'; else, visible = 'off'; end   
            set( hndls(j),'Color',[0 0 0],'Visible',visible);
        else
            % otherwise use the color mode
            switch( data.colormode )
             case 1, c = [ 0 0 0];
             case 2, c = data.spikes.info.kmeans.colors(clus,:);
             case 3,
                 id = data.spikes.assigns( find( data.spikes.info.kmeans.assigns == clus, 1 ) );
                 c = data.spikes.info.kmeans.colors(id,:);
            end
            set( hndls(j),'Color',c);
        end
    end
    
    % update legend
    if length(hndls) > 1 & ~data.show_density
      set_legend(data.spikes, data.colormode )
      legend hide
      data.show_legend = 0;
      set(hax,'UserData',data);
      set( findobj(get(hax,'UIContextMenu'),'Tag','legendoption'), 'Checked', 'off');
    end
end

% calculate a particular feature
function x = get_feature( spikes, which, feature, param )
  
% determine whether we are calculating a feature for outliers
is_outliers = which == 0;
if is_outliers
    data_source = spikes.info.outliers;
    which = 1:length(data_source.spiketimes);
else
    data_source = spikes;
end
    
switch (feature), 
    case 'Signal', 
        x = data_source.waveforms(which,param);      
    case 'PC', % (for principal components, we might need to compute the first time)
        if is_outliers
            x = data_source.waveforms(which,:) * spikes.info.pca.v(:,param);
        else
            x = spikes.waveforms(which,:) * spikes.info.pca.v(:,param);
        end
 
	case 'Cluster',
        if is_outliers, x = zeros( [1 length(which)] );
        else,           x = spikes.assigns(which); end
	case 'Minicluster',
        if is_outliers, x = zeros( [1 length(which)] );
        else,           x = spikes.info.kmeans.assigns(which); end		        
    case 'Event Time',
        x  = data_source.unwrapped_times(which); 
    case 'ISI Preceding'
        x = data_source.isis(which);
    case 'Total Energy',
        x = sum(data_source.waveforms(which,:).^2, 2).^.5;
        
    case 'Amplitude',  % again, we might need to compute these the first time they're used

        x = range( squeeze( data_source.waveforms(which,:,param) ),2 );
    case 'Width' % have to computer this as well
            w = data_source.waveforms(which,:,:);
            dc = spikes.info.detect.event_channel(which);
            [junk, posmin] = min(w,[],2);
            [junk, posmax] = max(w,[],2);
            widths = squeeze(abs(posmin-posmax));
            for j = 1:size(w,1)
               x(j) = widths(j, dc(j));       
            end
end
end

% callback to delete feature controls
function delete_controls(varargin)
     hax = findobj(gcf,'Tag','feature_plot');
     ssg = get(hax,'UserData');
     if isstruct(ssg)
      delete([ssg.xcontrol, ssg.ycontrol]);
     end
end
    
% callback to create a feature controller
function make_control(hObject, event, controlaxis)

    hax = findobj(gcf,'Tag','feature_plot');
    [o,h] = gcbo;

    data = get(hax,'UserData');
   
    controlhandle = data.([controlaxis 'control']);
    if (~isempty(controlhandle))
        figure(controlhandle);
    else
        data.([controlaxis 'control']) = feature_popup(h, hax, controlaxis);
        set( hax,'UserData',data)
    end
    uistack([data.xcontrol data.ycontrol], 'top');
end

% update the legend based on color mode
function set_legend( spikes, mode )

    % get axes for each cluster
    hndls = findobj(gca,'Type','line','Visible','on');
    hlist = []; clist = [];
    for j = 1:length(hndls)
        c = get(hndls(j),'UserData');
        if ~isempty(c), hlist(end+1) = hndls(j);  clist(end+1) = c; end 
    end
    [clist, i ] = sort(clist);
    cstrlist = {};  for j = 1:length(clist), cstrlist{j} = num2str(clist(j) ); end
    hlist(i) = hlist;
    
    if clist(1) == 0
        cstrlist{1} = 'outliers';
    end
    
    
    if mode == 1
        legend(hndls(1),'Data point')
        set( get(legend,'title'), 'String', 'All data')
    elseif mode == 2
        legend(hlist,cstrlist)
        set(get(legend,'title'),'String','Miniclusters')
    elseif mode == 3
       indices = find( ismember( clist, [unique(spikes.assigns(:))' 0] ) );
        legend(hlist(indices),cstrlist{indices})
        set(get(legend,'title'),'String','Clusters')       
    end
end
    
     
function spikes = add_isis_field( spikes )

    if isfield( spikes, 'assigns'), assigns = spikes.assigns;
    elseif isfield( spikes.info,'kmeans') assigns = spikes.info.kmeans.assigns;
    else  assigns = ones( size( spikes.unwrapped_times ) ); end
    spikes.isis = zeros( size(assigns) );
    
    % for normal spikes
    for j = unique( assigns )
        % select the spikes
        which = find( assigns == j );
        % calculate dt
        dt = diff( sort( spikes.unwrapped_times(which) ));
        dt = [mean(dt) dt];
        % save to data.isis
        spikes.isis(which) = dt;
    end
    maxdt = spikes.params.display.max_isi_to_display;
    spikes.isis( spikes.isis > maxdt )  = maxdt;

    % for outlier spikes
    if isfield( spikes.info, 'outliers' )
        dt = diff( sort(spikes.info.outliers.unwrapped_times ));
        spikes.info.outliers.isis = [mean(dt) dt];
        maxdt = spikes.params.display.max_isi_to_display;
        spikes.info.outliers.isis( spikes.info.outliers.isis > maxdt )  = maxdt;
    end
       
end
%%
%
%  FEATURE SELECT POP-UP CODE
%

% builds a feature controller
function newh = feature_popup(h, ax, which )
 
    % get mouse coordinates
    coords = get(h,'Units'); set(h,'Units','pixel');    fig_pos = get(h,'Position');ptr_pos = get(h,'CurrentPoint'); set(h,'Units',coords);
    
    % create figure at mouse coordinates
    newh = open('featureselect.fig');
    pos = get(newh,'Position');
    set(newh,'Position', [ fig_pos(1:2)+ptr_pos-[150 50] pos(3:4)],'DeleteFcn', {@popup_deleted, which, ax},'Name',[ which '-axis control'],'UserData',ax)
    
    % set values of controls
    data = get( findobj(h,'Tag','feature_plot'),'UserData');
    val = data.([which 'param']);
    choice = data.([which 'choice']);

    % popup
    list_box = findobj(newh,'Tag','popup_datatype');
    entries = get_popup_entries(data.spikes);
    for j = 1:length(entries), if isequal(entries{j},choice), target = j; end, end
    set( list_box,'String',entries,'Value',target,'Callback',{@listbox_callback,which});
    
    % value box
    param_box = findobj(newh,'Tag','edit_param1');    
    if isempty( val ), set( param_box,'Enable','off')
    else, set( param_box, 'String',num2str(val) ), end
    set(param_box,'Callback',{@param_callback,which});
    
    % advances
    set( findobj(newh,'Tag','button_last'),'Callback', {@control_click, -1, which});
    set( findobj(newh,'Tag','button_next'),'Callback', {@control_click, +1,which});
end

% callback when a feature is selected
function listbox_callback(hobject,event,which)
    
    [lb,h] = gcbo;
    
    ax = get(h,'UserData');
    data = get(ax,'UserData');
    entries = get(lb,'String');
    choice =  entries{ get(lb,'Value') };
    data.([which 'choice']) = choice;
  
    param_box = findobj(h,'Tag','edit_param1');     
    if which=='y' & isequal(choice,'PC'), defaultval = 2;  else, defaultval = 1; end
    if isequal(choice,'PC') | isequal( choice,'Signal') | isequal(choice,'Amplitude')
        set(param_box,'Enable','on','String',num2str(defaultval));
        data.([which 'param']) = defaultval;
    else
        set(param_box,'Enable','off','String','');
        data.([which 'param']) = [];
    end
    
    set(ax,'UserData',data);   
    redraw_features( ax );
end

% callback when a parameter value is updated
function param_callback(hobject,event,which)

    [a,h] = gcbo;
    list_box = findobj(h,'Tag','popup_datatype');
    entries = get(list_box,'String');
    datatype = entries( get(list_box,'Value') );
    
    param_box = findobj(h,'Tag','edit_param1');
    val = str2num( get(param_box,'String') );
    
    [lb,h] = gcbo;
    ax = get(h,'UserData');
    data = get( ax, 'UserData');
    
    % validate the value
    valid = rem(val,1) == 0 & val > 0;
    if isequal(datatype,'PC')
        valid = valid &  val <= size(data.spikes.waveforms(:,:), 2);
    elseif isequal( datatype,'Signal') 
        valid = valid &  val <= size(data.spikes.waveforms(:,:), 2);
    elseif isequal( datatype,'Amplitude')
        valid = valid & val <= size(data.spikes.waveforms,3);
    end
    
    if valid
          data.([which 'param']) = val;
          set( ax, 'UserData', data );
         redraw_features( ax);
    else
        beep
        set( param_box, 'String', num2str( data.([which 'param'] ) ) );
    end
end

% callback for when increment/decrement buttons are clicked for the feature control
function control_click(hObject,event,val,which)
   
    [a,h] = gcbo;
    param_box = findobj(h,'Tag','edit_param1');
    if isequal( get(param_box,'Enable'), 'on' )
      value = str2num( get(param_box,'String') );
      set(param_box,'String',num2str(value+val))
      param_callback([],[],which);
    end
end

% generates list of available features
function list = get_popup_entries( spikes)

    list = {'Signal','PC', 'Cluster', 'Minicluster', 'Event Time', 'ISI Preceding','Total Energy','Amplitude','Width'};
    strike = [];
    if ~isfield(spikes,'assigns'), strike = [strike 3]; end
    if ~isfield(spikes.info,'kmeans'), strike = [strike 4]; end
    list(strike) = [];
end

% tells plot that the control no longer exists
function popup_deleted(hObject,event, which, ax )
    data = get(ax,'UserData');
    data.([which 'control']) = [];
    set(ax,'UserData',data);
end    
   