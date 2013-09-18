function plot_waveforms( spikes, show, colormode,issubtree)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_waveforms - display waveforms for a collection of spikes
%
% Usage:
%           plot_waveforms( spikes, show, colormode,issubtree)
%
% Description:  
%   This function displays a collection of spike waveforms in a variety
% of ways.  There is both a colormode and a displaymode. Not all colormodes 
% and displaymodes are compatible and may not be available depending on how 
% the data has been processed.
%
%  Colormodes:  1 - all waveforms are displayed the same color
%               2 - all waveforms are displayed in different random colors
%               3 - all waveforms are displayed according to minicluster color
%               4 - all waveforms are displayed according to cluster color
%
%  Displaymodes: 1 - plot line for each waveform
%                2 - plot transparent error band for each minicluster or cluster (depends on colormode)
%                3 - plot 2d-histogram as a heat map (map set by spikes.params.display.cmap)
%
%  The default displaymode is set by spikes.params.display.default_waveformmode.  
%  The current colormode and displaymode can be changed by right clicking the axes 
%  and selcting from a context menu.  There is also the option to impose the 
%  current set of mods on all other waveform plots in the same figure.  Some
%  modes also allow for the display of a legend that can be accessed from
%  the same menu.
%
%  When using display mode 1 or 2 with color mode 3 or 4, left clicking on 
%  data will raise all data of the same color to the front.  Right clicking 
%  will send it to the back.
%
%  A time scalebar is included in the lower right hand corner of the plot.
%  Its duration is set by spikes.params.display.time_scalebar.  Its color
%  is made to match the cluster or minicluster color if possible.
%
%  Multi-channel waveforms are concatenated and red vertical lines are used to
%  separate the waveform on different channels.
%
%  The y-axis attempts to identify the cluster or minicluster being plotted
% unless the data is mixed from different miniclusters.  Also specified are
% the total number of spikes and the number of refractory period vioaltions
% (RPVs) in the unwrapped spike time as defined by 
% spikes.params.refractory_period.
%
% Input:
%   spikes - a spike structure
%
% Optional input:
%   show          - array describing which events to show in plot
%                 - see get_spike_indices.m, (default = 'all')
%   colormode     - mode for coloring individual data points 
%                 - 1=> all different, 2=>all same, 3=>use minicluster color, 4=>use cluster color
%                 - the default is to use the highest flag that is valid for the given spikes object
%   issubtree     - mainly for use by split_tool.m to alter the behavior of 
%                   this function when specifically displaying a sub branch
%                   of the full aggregration tree that is not itself a cluster
%                  - default = 0
   
    cla reset;
    % argument checking
   if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end
   if nargin< 4, issubtree = 0; end
   valid_modes = [1 isfield(spikes.info,'kmeans') isfield(spikes,'assigns')];
   if nargin < 3, colormode =  1 + find(valid_modes,1,'last' ); end
   if nargin < 2, show = 1:size(spikes.waveforms,1); end
    
    % which spikes are we showing?
    show = get_spike_indices(spikes, show );

    % populate data structure 
    data.ylims = [min(spikes.waveforms(:)) max(spikes.waveforms(:))];
    data.valid_modes = valid_modes;
    data.colormode = colormode;   
    spiketimes =  sort( spikes.unwrapped_times(show) );
    data.rpv  = sum( diff(spiketimes)  <= (spikes.params.refractory_period * .001) );
    data.waveforms = spikes.waveforms(show,:,:);
    data.cmap = spikes.params.display.cmap;
    if any(data.valid_modes(2:3))
        data.colors = spikes.info.kmeans.colors;
        if data.valid_modes(3), data.assigns   = spikes.assigns(show); end
        if data.valid_modes(2), data.subassigns   = spikes.info.kmeans.assigns(show); end
    end
    data.time_scalebar =  spikes.params.display.time_scalebar;
    data.displaymode = spikes.params.display.default_waveformmode;
    data.defaultcolor = [0.5 0 0];
    data.Fs = spikes.params.Fs;

    % set up title and bar color
    
    % 1 cluster and coloring by cluster
    if data.colormode == 4 & length(unique(spikes.assigns(show))) == 1
        data.title = ['Cluster # ' num2str( spikes.assigns(show(1)) ) ];
        data.barcolor = data.colors(  spikes.assigns(show(1)),: );
    
    % this is a bit of a hack:  if we are plotting a subtree, want to make color
    % and name same as subtree tope node.         
    elseif issubtree
        subclusterlist = unique( spikes.info.kmeans.assigns(show) )';
        idx = find( ismember( spikes.info.tree(:,1), subclusterlist ), 1,'last');
        if isempty(idx)
            subclus =  spikes.info.kmeans.assigns(show(1));
        else
            subclus = spikes.info.tree(idx,1);
        end
        data.title = ['Minicluster # ' num2str(subclus) ];
        data.barcolor = data.colors( subclus,: );
     
    % 1 minicluster and coloring by minicluster
    elseif data.colormode == 3 & length(unique(spikes.info.kmeans.assigns(show))) == 1
        data.title = ['Minicluster # ' num2str( spikes.info.kmeans.assigns(show(1)) ) ];
        data.barcolor = data.colors( spikes.info.kmeans.assigns(show(1)),: );
        
    % multiple clusters
    elseif data.valid_modes(3) 
        data.title = [num2str( length(unique(spikes.assigns(show))) ) ' clusters'];
        data.barcolor = [1 0 0];
    
    % multiple miniclusters
    elseif  data.valid_modes(2) 
        data.title = [num2str( length(unique(spikes.info.kmeans.assigns(show))) ) ' miniclusters'];
        data.barcolor = [1 0 0];
        
    % unclustered data
    else
        data.title = 'Unclustered data';
        data.barcolor = [1 0 0];
    end
        
    data.subcluster = data.colormode - 3;
    
    % save user data
    set(gca,'UserData', data,'Tag','waveforms' )

    % update display
    update_waveforms( [], [], data.colormode, data.displaymode, gca);

% callback for updating appearance of waveforms
function update_waveforms( hObject, event, colormode, displaymode, ax)
   
    if displaymode == 2, colormode = max(colormode,2); end
    if displaymode == 3, colormode = 2; end
    
    data = get( ax,'UserData');
    data.colormode = colormode;
    data.displaymode = displaymode;
    set(ax,'UserData',data)
    set(gcf,'CurrentAxes',ax )
    
    % call plot function
    make_waveforms(ax, data);

    set(ax,'UIContextMenu', get_menu(colormode,displaymode,data.valid_modes, ax) )

% callback to switch color and display mode for all axes
function impose_all( hObject, event, colormode, displaymode, ax)
   
       [o,h] = gcbo;     
   
        my_axes = findobj( h,'Tag','waveforms');
        my_axes = setdiff( my_axes, ax );
        for j = 1:length(my_axes)
                update_waveforms( [], [], colormode, displaymode, my_axes(j));
        end

% make user context menu    
function cmenu = get_menu( colormode, displaymode,valid_modes,ax )

    cmenu = uicontextmenu;

    % COLORMODES
    colormodes = {'All different','All same','By minicluster', 'By cluster'};
    if ~valid_modes(3), colormodes(4) = []; end
    if ~valid_modes(2), colormodes(3) = []; end
        
    for j = 1:length(colormodes)
        c(j) = uimenu(cmenu, 'Label', colormodes{j}, 'Callback', {@update_waveforms, j, displaymode,ax} );
    end
    set(c(colormode),'Checked','on');
    if displaymode == 2, set( c(1),'Enable','off'); end
    if displaymode == 3, set( c,'Enable','off'); end
    
    %DISPLAYMODES
    d(1) = uimenu(cmenu, 'Label', 'Raw', 'Callback', {@update_waveforms, colormode, 1,ax},'Separator','on');
    d(2) = uimenu(cmenu, 'Label', 'Bands', 'Callback',{@update_waveforms, colormode, 2,ax});
    d(3) = uimenu(cmenu, 'Label', '2D histogram', 'Callback',{@update_waveforms, colormode, 3,ax});
    set(d(displaymode),'Checked','on');

    %IMPOSE ON ALL
    uimenu(cmenu, 'Label', 'Use this style on all waveforms in figure', 'Callback', {@impose_all, colormode, displaymode,ax},'Separator','on');

    %LEGEND
    if colormode >= 3 & displaymode < 3
      uimenu(cmenu, 'Label', 'Show legend', 'Callback',@toggle_legend,'Separator','on','Checked',get(legend,'Visible'),'Tag','legend_option');
    end
    
% callback to hide/show legend
function toggle_legend(varargin)
     
    item = findobj( get(gca,'UIContextMenu'),'Tag', 'legend_option');
    
    if isequal( get(legend,'Visible'),'on')
       legend('hide'); set(item,'Checked','off');
   else
       legend('show'); set(item,'Checked','on');
    end
   
% plot to didsplay waveforms
function make_waveforms(ax,data)

% unpack
displaymode = data.displaymode;
colormode = data.colormode;
waveforms = data.waveforms;
defaultcolor = data.defaultcolor;
time_scalebar = data.time_scalebar;
Fs           = data.Fs;

% interpet colormodes
if data.colormode == 2
    colors = defaultcolor;
    assigns = ones([1 size(waveforms,1)]);
 elseif data.colormode == 3
    colors = data.colors;
    assigns = data.subassigns;
elseif data.colormode == 4
    colors = data.colors;
    assigns = data.assigns;
end

% prep axes
cla; 
legend off
set(gca,'Color',[ 1 1 1])
hold on

%
% DISPLAY WAVEFORMS
%
num_samples = size(waveforms(:,:),2);

% plot all waveforms
if displaymode == 1 & colormode == 1
       plot( 1:size(waveforms(:,:),2), waveforms(:,:) );
       
% use 2D density plot
elseif displaymode == 3
     cmap = data.cmap;
     [n,x,y] = histxt(waveforms(:,:));  
     h = imagesc(x,y,n);  
     colormap(cmap);
     set( gca,'Color', cmap(1,:) );    

% plot waveforms in groups
else
    clusts = sort(unique(assigns));
    for j = 1:length(clusts)
      
        color = colors( clusts(j),:);
        mine = find( assigns == clusts(j) );
        
        % plot lines
        if displaymode == 1
            lh(j) = mplot(1:num_samples, waveforms(mine,:), 'Color', color);
            set(lh(j), 'ButtonDownFcn', {@raise_me});

        % plot standard deviation bands
        elseif displaymode == 2
            
             % group traces and show +/- 2 standard deviations            
              [lh(j) ,ph(j)] = error_area(mean(waveforms(mine,:),1), 2*std(waveforms(mine,:),1,1));
              set(lh(j), 'Color', brighten(color, -0.6), 'ZData', clusts(j)* ones( size(get(lh(j),'XData')))  );
              set(ph(j), 'FaceColor', color, 'ZData', clusts(j)* ones( size(get(ph(j),'XData'))), 'FaceAlpha', 0.8);  
              set([lh(j) ph(j)], 'ButtonDownFcn', {@raise_band, [lh(j) ph(j)] });        
           
        else
              error(['Invalid plot_waveforms display mode (' num2str(displaymode) ') or color mode (' num2str(colormode) ').' ] );
        end
    end
end
hold off

set(gca,'Tag','waveforms','XLim',[1 num_samples],'YLim',data.ylims)

% make vertical lines to separate wavform into its channels
num_channels = size(waveforms,3);
num_samples = size(waveforms,2);
ylims = get(gca,'YLim');
if num_channels > 1
    if displaymode == 2
       for j = 1:num_channels-1
            l(j) = line( 1 + num_samples * j * [1 1], ylims, max(clusts)*[1 1] + 1);           
       end
        set(l,'ButtonDownFcn', {@raise_band, l})
    else
       for j = 1:num_channels-1
            l(j) = line( 1 + num_samples * j * [1 1], ylims);
       end
           set(l,'ButtonDownFcn', {@raise_me})
    end
      set( l, 'Color',data.barcolor,'LineWidth',1.5 ) % electrode dividers
end

% scale bar
ms =  time_scalebar*Fs/1000;
maxX = num_channels*num_samples;
lineY = ylims(1) + (ylims(2)-ylims(1)) * .075;
l = line(  maxX - 5 - [0 ms], lineY*[1 1] );

% set the raise callback
if displaymode == 2,
    set(l,'ButtonDownFcn', {@raise_band, l})
else
       set(l,'ButtonDownFcn', {@raise_me})
end
set(l,'Color',data.barcolor,'LineWidth', 3)

%  make legend only if multiple colors are used
if colormode >= 3 & displaymode < 3 
    leg = cell(length(clusts),1);
   for k = 1:length(clusts),  leg{k} = num2str(clusts(k));  end;
   l = legend(lh,leg,'Location','Best');
   set(l,'FontSize',7);
end
legend hide

% label axes
xlabel('Sample')
ystr = { data.title, ['N = ' num2str(size(waveforms,1)) '  (' num2str(data.rpv) ' RPVs)'] };
ylabel(ystr)
