function splitmerge_tool( spikes, show_clusters, h, varname)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% splitmerge_tool - tool for inspection and correction of a sorted spikes object
%
% Usage:
%      splitmerge_tool( spikes, show_clusters, h, varname)
%
% Description:  
%   This is a flexible tool for viewing clusters, merging or splitting them,
% labeling their status, and removing outliers. See spike_sorting_manual.pdf 
% for instructions in how to use this tool.
%   
%  See also the other tools which splitmerge_tool uses:
%          split_tool
%          outlier_tool
%          sliderFigure
%
% See also ss_default_params.m as several fields in spikes.params.display affect
% the appearance of this tool.
%
%             default_figure_size  - size of figure
%             figure_font_size    - font size
%             merge_fig_color      - color of figure background
%             label_categories     - labels that can be applied to clusters
%             label_colors         - corresponding colors for the panels of different labels
%             
% Input: 
%   spikes - a spikes object
%
% Optional Input:
%   show_clusters  - list of IDs for clusters to be worked with (default is to use all clusters)
%   h              - handle to figure to be used for the splitmerge tool (default = gcf)
%   varname        - default variable name to which to save modified spikes object (default = 'spikes')
%
   
    %  check arguments
    if ~isfield(spikes,'assigns'), error('No assignments found in spikes object.'); end
    if nargin < 3, h = figure('Units','Normalized','Position',spikes.params.display.default_figure_size); end
    if nargin < 2 | isequal( show_clusters,'all') 
        figdata.hidden_clusters = []; 
    else
        figdata.hidden_clusters = setdiff( unique(spikes.assigns), show_clusters );
    end
    figdata.master_hidden = figdata.hidden_clusters;
    if nargin < 4, varname = 'spikes'; end
    
    % set up figure
    clf(h,'reset');
    set(h,'defaultaxesfontsize',spikes.params.display.figure_font_size);    
    set(h,'Interruptible','off','BusyAction','cancel')
    set(h,'Color',spikes.params.display.merge_fig_color);
    set(h,'KeyPressFcn',@hot_keys)  % hot key callback

    % make toolbar buttons
    th = uitoolbar(h);
    data.ims = load('icons_ims');
    % saving and loading
    uipushtool(th,'CData',data.ims.im_save,'ClickedCallBack',{@saveSpikes,h},'userdata',varname,'Tag','saveButton','TooltipString','Save to workspace (s)','Separator','on');
    figdata.sfb = uipushtool(th,'CData',data.ims.im_savefile,'ClickedCallBack',{@saveSpikesToFile,h},'Tag','saveFileButton','TooltipString','Save to file');
    uipushtool(th,'CData',data.ims.im_loadfile,'ClickedCallBack',{@loadSpikesFromFile,h},'Tag','loadFileButton','TooltipString','Load from file');
    
    % calling the other figures
    uipushtool(th,'CData',data.ims.im_show,'ClickedCallBack',{@makeShowClust,h},'TooltipString','Show selected clusters in separate figure','Separator','on');
    uipushtool(th,'CData',data.ims.im_sep,'ClickedCallBack',{@makeSeparationAnalysis,h},'TooltipString','Compare selected cluster in separate figure');
    uipushtool(th,'CData',data.ims.im_scatter,'ClickedCallBack',{@makePlotFeatures,h},'TooltipString','Plot features of selected clusters in separate figure');
    
    % show and execute merge
    uipushtool(th,'CData',data.ims.im_eye,'ClickedCallBack',{@show_selected,h},'TooltipString','Show plots for merger of selected clusters (e)','Separator','on');
    uipushtool(th,'CData',data.ims.im_go,'ClickedCallBack',{@execute_change,h},'Tag','executeButton','TooltipString','Merge selected clusters (x)');

    % manipulation of panels
    uipushtool(th,'CData',data.ims.im_select,'ClickedCallBack',{@select_all,h},'TooltipString','Select all panels (a)','Separator','on');
    uipushtool(th,'CData',data.ims.im_deselect,'ClickedCallBack',{@deselect_all,h},'TooltipString','Deselect all panels (d)');
    uipushtool(th,'CData',data.ims.im_hide,'ClickedCallBack',{@hide_selected,h},'TooltipString','Hide selected panels (h)');
    uipushtool(th,'CData',data.ims.im_reveal,'ClickedCallBack',{@reveal_all,h},'TooltipString','Reveal hidden panels (r)','Enable','off','Tag','reveal_button');
    uipushtool(th,'CData',data.ims.im_rearrange,'ClickedCallBack',{@rearrange_panels,h},'TooltipString','Match panel tiling to figure size (p)','Enable','on');

    % save spikes data to figure
    figdata.spikes = spikes;
    figdata.selected  = [];
    figdata.clus_list = sort(unique(spikes.assigns));
    figdata.handles.panels = []; figdata.handles.axes = []; figdata.handles.id=[];
    
    % initialize the panels
    p1 = make_uipanel_on_grid( [1 1 4 1], num2str(['Combined selected clusters']), spikes.params.display, h );
    figdata.handles.cluster_panel = p1;
    set(h,'UserData',figdata)
    init_panels(h);
    
    % initialize the slider figure
    sliderFigure( h, spikes.params.display.outer_margin, th, spikes.params.display.aspect_ratio )

    % set up figure
    set(h,'Name','Merge Tool','NumberTitle','off')
    figure(h)
    
end

%%
%  callbacks for toolbar buttons
%

% save spikes object to workspace
function saveSpikes(varargin)
    
    h = varargin{end};
    figdata = get(h,'UserData');
    me = findobj(h,'Tag','saveButton');
    
    % ask for name of output variable
    varname = inputdlg('Save spike data in what variable?', 'Save to workspace', 1,{get(me,'UserData')});
    if ~isempty(varname)
        assignin('base',varname{1},figdata.spikes)
        set(me,'UserData',varname{1} );
    end
end

% save spikes object to file
function saveSpikesToFile(varargin)

    h = varargin{end};
    figdata= get(h,'UserData');
    b = figdata.sfb;
    
    % check if a filename was previously selected
    filename_default = get(b,'UserData');
    if isempty(filename_default)
        filename = 'spikes.mat';
    else
        filename = [filename_default.pathname filename_default.filename];
    end
    
    % open up dialog to ask for file location
    [FileName,PathName,FilterIndex] = uiputfile('*.mat','Save Spikes object',filename);
    if ~isequal(FileName,0)
        
        spikes = figdata.spikes;
        save([PathName FileName],'spikes');

        % save filename as new 
        a.pathname = PathName;
        a.filename = FileName;
        set(b,'UserData',a);
    end
end
     
% load spikes object from file
function loadSpikesFromFile(varargin)
    
    h = varargin{end};
      % look for a previously used file name
    filename_default = get(findobj(h,'Tag','saveFileButton'),'UserData');
    if isempty(filename_default)
        filename = [pwd '\spikes.mat'];
    else
        filename = [filename_default.pathname filename_default.filename];
    end
    [filename,pathname] = uigetfile('*.mat','Load spikes object from file.',filename);
    
    if ~isequal(filename,0)
      
      % load file
      a = load([pathname filename]);  
      names = fieldnames(a);
      spikes = getfield(a,names{1});
      
      % inefficient --> hard reset the tool
      clf(h,'reset')
      splitmerge_tool(spikes, unique(spikes.assigns),h)
     
    end
end
     
% make instance of show_clusters
function makeShowClust(varargin)

      figdata = get(varargin{end},'UserData');
      selected = figdata.selected;
      
      if isempty(selected)
                show_clusters(figdata.spikes);
      else
                show_clusters(figdata.spikes,selected);
      end
end

% make instance of compare_clusters
function makeSeparationAnalysis(varargin)
      figdata = get(varargin{end},'UserData');
      selected = figdata.selected;
      if isempty(selected)
          compare_clusters(figdata.spikes);
      else
          compare_clusters(figdata.spikes,selected);
      end
end

% make instance of plot_features
function makePlotFeatures(varargin)
      figdata = get(varargin{end},'UserData');
      figure
      selected = figdata.selected;
      if isempty(selected), selected = 'all'; end
      plot_features(figdata.spikes,selected)
end

% show plots for merge of selected panels
function show_selected(varargin)

    h = varargin{end};
    figdata = get(h,'UserData');
    clear_display_panel( h )
    
    if ~isempty( figdata.selected )
        display_cluster( h, figdata.selected, 1, figdata.handles.cluster_panel, 1 );
    end

end

% execute a merge on the selected panels
function execute_change(varargin)

    % get figure data
    h = varargin{end};
    figdata = get(h,'UserData');
      
    % must select at least two clusters
    if length(figdata.selected) < 2
        warndlg( 'Select two or more clusters to merge', 'Not enough clusters selected');
        
    else        
        % get list of all subclusters to remove
        for j = 1:(length(figdata.selected)-1)
            figdata.spikes = merge_clusters( figdata.spikes, figdata.selected(end), figdata.selected(j) );
        end
        figdata.clus_list = sort(unique(figdata.spikes.assigns));
        kept_clus =  intersect(figdata.selected, figdata.clus_list );
        set(h,'UserData',figdata);
      
        % out with the old
        delete( get(figdata.handles.cluster_panel,'Children') )
        gone = setdiff(figdata.selected, kept_clus );
        which = find( ismember( figdata.handles.id, gone ) );
        delete( figdata.handles.panels(which) );
        figdata.handles.panels(which) = [];
        figdata.handles.axes(which) = [];
        figdata.handles.id(which) = [];
       
        % fix the new  
        which = find( kept_clus == figdata.handles.id );
        delete( figdata.handles.axes(which));
        a = axes('Parent',figdata.handles.panels(which));
        set(h,'CurrentAxes',a)
        plot_waveforms( figdata.spikes, figdata.handles.id(which) );
        figdata.handles.axes(which) = a;  

        % reset panel color and position
        color = figdata.spikes.params.display.label_colors(1,:);
        set(figdata.handles.panels(which),'BackgroundColor',color);        
        
        figdata.selected = [];
        set(h,'UserData',figdata)
        rearrange_panels(h);
    end
    
end          

% select all panels
function select_all( varargin )
    h = varargin{end};
    figdata = get(h,'UserData');
    set(figdata.handles.panels,'BackgroundColor',[1 1 1] );
    figdata.selected = setdiff(figdata.handles.id,figdata.hidden_clusters);
    set(h,'UserData',figdata)   
end

% deselect all panels
function deselect_all( varargin )
    h = varargin{end};
    figdata = get(h,'UserData');
    for j = 1:length(figdata.handles.panels)
        p = figdata.handles.panels(j); 
        clus = str2num(get(p,'Title'));
        idx = figdata.spikes.labels(figdata.spikes.labels(:,1)==clus,2);
        color = figdata.spikes.params.display.label_colors(idx,:);
        set(p,'BackgroundColor',color);
    end
    figdata.selected = [];
    set(h,'UserData',figdata)   
end

% hide selected panels
function hide_selected(varargin)  

    h = varargin{end};
    figdata = get(h,'UserData');
    if ~isempty( figdata.selected)
      figdata.hidden_clusters = [figdata.hidden_clusters figdata.selected ];
      
      % recolor the selected
      for j = 1:length(figdata.selected)
        which1 = find(figdata.handles.id==figdata.selected(j) );
        which2 = find(figdata.spikes.labels(:,1) ==figdata.selected(j));
        color = figdata.spikes.params.display.label_colors(figdata.spikes.labels(which2,2),:);
        set(figdata.handles.panels(which1),'BackgroundColor',color);
      end
      figdata.selected = [];
       
       set( findobj(h,'Tag','reveal_button'), 'Enable','on');
      set(h,'UserData',figdata) 
      rearrange_panels( h );
    end
end
 
% rearrange panels (usually called after resizing figure)
function rearrange_panels( varargin )

    h = varargin{end};
    % make all panels and axes invisible
    figdata = get(h,'UserData');
    set( figdata.handles.panels, 'Visible','off')
    
    % get list of clusters we are supposed to show
    show_list = sort( setdiff( unique(figdata.spikes.assigns), [figdata.hidden_clusters figdata.master_hidden] ) );
    
    % move them to correct location
    fpos = get(h,'Position');
    d= figdata.spikes.params.display;
    r_button = findobj(h,'Tag','reset');
    if isempty(r_button), zoom_factor = 1; else, zoom_factor = get(r_button,'UserData'); end
    num_cols = max(1,floor( (fpos(3) - d.outer_margin) / ((d.margin+d.width)/zoom_factor) ));
    for j = 1:length(show_list)
        % where would this one go?
        row = ceil(j/num_cols);
        col = rem(j,num_cols);  if col==0,col=num_cols;end
    
       % find which one
       a = find( figdata.handles.id == show_list(j) ); 
       
       % set new position
       set(figdata.handles.panels(a), 'Position', get_pos_on_grid([1+row col 1 1], figdata.spikes.params.display, h ) );
       
    end
    % make them visible
    indices = find( ismember( figdata.handles.id, show_list ) );
    set( figdata.handles.panels(indices),'Visible','on')
    
end

% reveal panels
function reveal_all(varargin)
    h = varargin{end};
    figdata = get(h,'UserData');
    figdata.hidden_clusters = figdata.master_hidden;
    set(h,'UserData',figdata); 
    rearrange_panels( h );
    set( findobj(h,'Tag','reveal_button'), 'Enable','off');

end

% hot key definitions for quick access to these callbacks 
function hot_keys(varargin)

    h =varargin{end-1};
    event = varargin{end};
    switch(event.Key)
        case {'a'}, select_all(h)
        case {'d'}, deselect_all(h)
        case {'s'}, saveSpikes(h)
        case {'x'}, execute_change(h)
        case {'e'}, show_selected(h)
        case {'h'}, hide_selected(h)
        case {'r'}, reveal_all(h)
        case {'p'}, rearrange_panels(h);

    end
end

%% 
%  Panel functions  
%

% initialize all panels
function init_panels( h )

    figdata = get(h,'UserData');
     set(h,'Pointer','watch')
     pause(.01)
    figdata.clus_list = setdiff( sort(unique(figdata.spikes.assigns)), figdata.hidden_clusters );
     
    for k = 1:length(figdata.clus_list)
        
        % create waveform panel
        p = uipanel('Title',num2str(figdata.clus_list(k)),'Units','pixels','Visible','off');
        set(p,'Title',num2str(figdata.clus_list(k)) );
        set(p,'ButtonDownFcn',{@click_panel,h, p} );
        idx = figdata.spikes.labels( figdata.spikes.labels(:,1) == figdata.clus_list(k), 2 );
        color = figdata.spikes.params.display.label_colors(idx,:);
        if ismember( figdata.clus_list(k), figdata.selected), color = [1 1 1]; end
        set(p,'BackgroundColor',color);
        set(p,'UIContextMenu',get_panel_menu( figdata.spikes.params.display.label_categories, idx, p, h ) );  
        
        % create waveform axes
        a = axes('Parent',p);   
        plot_waveforms( figdata.spikes, figdata.clus_list(k) );
        
        % save info
        figdata.handles.panels(k) = p;
        figdata.handles.axes(k) = a;
        figdata.handles.id(k) = figdata.clus_list(k);
      
    end
    
    set( h, 'UserData',figdata);
    rearrange_panels(h);
    set(h,'Pointer','arrow')
   
end

% select/de-select a panel
function click_panel( varargin )
    h = varargin{end-1};
    p  = varargin{end};
    
    figdata = get(h,'UserData');
        clus = str2num(get(p,'Title'));
    
    % update list
    if ~isequal( get(h,'SelectionType'), 'alt' )
        if ismember( clus, figdata.selected )
            idx = figdata.spikes.labels(figdata.spikes.labels(:,1)==clus,2);
            color = figdata.spikes.params.display.label_colors(idx,:);
            set(p,'BackgroundColor',color);
            figdata.selected( figdata.selected == clus ) = [];
        else        
            set(p,'BackgroundColor',[1 1 1] );
            figdata.selected( end +1  ) = clus;
        end
        set(h,'UserData',figdata)
    end
end

% delete a set of panels
function delete_panels( ps )
    for j = 1:length(ps)
       d = get(ps(j),'UserData');
       delete( d.my_axes );
       delete(ps(j));
    end
end

% clear the main display panel
function clear_display_panel( h )
    
    figdata = get(h,'UserData');
    pdata = get(figdata.handles.cluster_panel,'UserData');
    delete( pdata.my_axes ); 
    pdata.my_axes = []; 
    set( figdata.handles.cluster_panel,'UserData',pdata);
    
end

%%
%  functions for panel context-menu
%

% create the menu
function cmenu = get_panel_menu( categories, me,p,h)
        cmenu = uicontextmenu;
        clus = str2num( get(p,'Title' ) );
        d(1) = uimenu(cmenu, 'Label', 'Open in split tool', 'Callback', {@open_split_tool,h,clus} );
        d(2) = uimenu(cmenu, 'Label', 'Open in outlier tool', 'Callback', {@open_outlier_tool,h,clus} );
        
        for j = 1:length(categories)
            d(j+2) = uimenu(cmenu, 'Label', categories{j}, 'Callback', {@change_label, p, j,h});
            if j == me, set( d(j+2), 'Checked','on','Enable','off'); end
        end
        set(d(3),'Separator','on');
end

% open instance of the split tool for the current panel
function open_split_tool( varargin )
    
    h = varargin{end-1};
    clus = varargin{end};
    figdata = get(h,'UserData');
    split_tool( figdata.spikes, clus,[],@remote_update,h );
end
    
% open instance of the outlier tool for the current panel
function open_outlier_tool( varargin )
    h=varargin{end-1};
    clus = varargin{end};
    figdata = get(h,'UserData');
    outlier_tool( figdata.spikes, clus, [], @remote_update,h)

end

% select a new label
function change_label(varargin)
    h = varargin{end};
    p = varargin{end-2};
    which = varargin{end-1};  
    clus = str2num(get(p,'Title'));
    figdata = get(h,'UserData');
    
    % change labels in spikes object
    idx = find( figdata.spikes.labels(:,1) == clus );
 
    figdata.spikes.labels(idx,2) = which;
    set(h,'Userdata',figdata);
    
    % change color of current cluster if not selected
    if ~ismember( clus, figdata.selected )
     set(p,'BackgroundColor', figdata.spikes.params.display.label_colors(which,:) );
    end   
    % update context menu
    set(p,'UIContextMenu',get_panel_menu(figdata.spikes.params.display.label_categories,which,p,h) );
    
end

% update this tool after spikes object is modified by split tool or the outlier tool
function remote_update( h, spikes, which )

    set(0,'CurrentFigure',h)
    figdata = get(h,'UserData');
    figdata.spikes = spikes;
    figdata.selected = [];
   
    %delete the upper panel
    delete( get( figdata.handles.cluster_panel,'Children'))
    
    % update the one we left for
    a = find( figdata.handles.id == which );
   
    if isempty(find(figdata.spikes.assigns==a))
        delete( figdata.handles.panels(a) );
        figdata.handles.panels(a) = [];
        figdata.handles.axes(a) = [];
        figdata.handles.id(a) = [];
    else
        delete( figdata.handles.axes(a));
        ax = axes('Parent',figdata.handles.panels(a)); 
        set(h,'CurrentAxes',ax)
        plot_waveforms( figdata.spikes, figdata.handles.id(a) );
        figdata.handles.axes(a) = ax;  
    end
 
    % get list of new clusters
    new_clusters = setdiff( unique(figdata.spikes.assigns), union( figdata.handles.id, figdata.hidden_clusters ) );
    pos = get_pos_on_grid([2 1 1 1], figdata.spikes.params.display, h ); 
    for j = 1:length(new_clusters)
        % create waveform panel
        p = uipanel('Title',num2str(new_clusters(j)),'Units','pixels','Visible','off','Position',pos);
        set(p,'Title',num2str(new_clusters(j)));
        set(p,'ButtonDownFcn',{@click_panel,h, p} );
        color = figdata.spikes.params.display.label_colors(1,:);
        set(p,'BackgroundColor',color);
        set(p,'UIContextMenu',get_panel_menu( figdata.spikes.params.display.label_categories, 1, p,h ) );  
        
        % create waveform axes
        a = axes('Parent',p);   
        plot_waveforms( figdata.spikes, new_clusters(j) );
           
        
        %build them
        figdata.handles.panels(end+1) = p;
        figdata.handles.axes(end+1) = a;
        figdata.handles.id(end+1) = new_clusters(j);        
    end
    
    %   sort those lists
    [figdata.handles.id,a] = sort( figdata.handles.id);
    figdata.handles.axes = figdata.handles.axes(a);
    figdata.handles.panels = figdata.handles.panels(a);
      
    
    % recolor panels
    for j = 1:length(figdata.handles.id)
        b = find( figdata.spikes.labels(:,1) == figdata.handles.id(j) );
        c = figdata.spikes.labels(b,2);
        color = figdata.spikes.params.display.label_colors(c,:);
        set(figdata.handles.panels(j),'BackgroundColor',color)
    end
    set(h,'UserData',figdata);

    % clear the display panel
    rearrange_panels(h);  
    figure(h)
end

%% 
%  utility functions
%

% populate the main display cluster with the merge of the selected panels
function display_cluster( h, clus, row, p, issubtree )
        
    % make 4 plots
    if nargin < 5, issubtree = 0; end
    
    figdata = get(h,'UserData');
    spikes = figdata.spikes;
    set(0,'CurrentFigure',h)
    
    % plot waveforms
    ax(1) = subplot(1,4,1,'Parent',figdata.handles.cluster_panel);
    if issubtree
        plot_waveforms(spikes,clus,3);
    else
        plot_waveforms(spikes,clus);
    end
    
    % plot detection criterion
    ax(2) = subplot(1,4,2,'Parent',figdata.handles.cluster_panel);
    plot_detection_criterion(spikes,clus);

    % plot inter-spike interval histogram
    ax(3) = subplot(1,4,3,'Parent',figdata.handles.cluster_panel);
    plot_isi(spikes,clus);
    
    % stability plot
    ax(4) = subplot(1,4,4,'Parent',figdata.handles.cluster_panel);
    [ax(4),ax(5)] = plot_stability(spikes,clus);
    
    % adjust margins -- so arbitrary!
    p = get(ax(1),'Position'); m = p(1)/3; set(ax(1),'Position',[p(1)-m p(2:4)]);
    p = get(ax(2),'Position'); set(ax(2),'Position',[p(1)-(m/3) p(2:4)]);
    p = get(ax(3),'Position'); set(ax(3),'Position',[p(1)+(m/3) p(2:4)]);
    p = get(ax(4),'Position'); set(ax([4 5]),'Position',[p(1)+m p(2:4)]);
    
end
