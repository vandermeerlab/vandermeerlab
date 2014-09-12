function split_tool( spikes, clus, h, savename, target_fig)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% split_tool - tool for inspection and division of a single cluster from a sorted spikes object
%
% Usage:
%      split_tool( spikes, clus, h, savename, target_fig)
%
% Description:  
%      Tool to inspect and break up a single cluster.  This tool can be called
% alone but is typically accessed through splitmerge_tool. See spike_sorting_manual.pdf 
% for instructions in how to use this tool.
%
% See ss_default_params.m as several fields in spikes.params.display affect
% the appearance of this tool.
%
%             default_figure_size          - size of figure
%             figure_font_size            - font size
%             split_fig_color              - color of figure background
%             label_colors                 - 1st entry sets color of panels
%             initial_split_figure_panels  - number of sub-trees to show by default
%             
% Also see sliderFigure.m since this tool makes use of sliderFigure controls.
%
%
% Input: 
%   spikes         - a spikes object
%   clus           - ID for the cluster to examine
%
% Optional Input:
%   h              - handle for figure to be used for the splitmerge tool (default = gcf)
%   savename       - default name for workspace variable to save to (default = 'spikes')
%                  - this can also be a function handle to call when saving
%   target_fig     - handle for splitmerge_tool figure that called split_tool (default = [])
%
   
    % check arguments
    if ~isfield(spikes,'assigns'), error('No assignments found in spikes object.'); end
    if nargin < 2, error('No cluster specified.'); end
    if nargin < 3, h = figure('Units','Normalized','Position',spikes.params.display.default_figure_size); end
    if nargin < 4, savename = 'spikes'; end
    if nargin < 5, target_fig = []; end
   
    % set up figure
    if isempty(h), h = figure('Units','Normalized','Position',spikes.params.display.default_figure_size); end
    clf(h,'reset');
    set(h,'defaultaxesfontsize',spikes.params.display.figure_font_size);
    set(h,'Interruptible','off','BusyAction','cancel')
    set(h,'Color',spikes.params.display.split_fig_color);
    set(h,'KeyPressFcn',@hot_keys)
    
    % make toolbar buttons
    th = uitoolbar(h);
    data.ims = load('icons_ims');
    
    uipushtool(th,'CData',data.ims.im_save_and_close,'ClickedCallBack',{@saveSpikes,h},'TooltipString','Save to merge tool and close (s)','Separator','on');

    % show selected panels using separate figure functions
    uipushtool(th,'CData',data.ims.im_show,'ClickedCallBack',{@makeShowClust,h},'TooltipString','Show selected clusters in separate figure','Separator','on');
    uipushtool(th,'CData',data.ims.im_sep,'ClickedCallBack',{@makeSeparationAnalysis,h},'TooltipString','Compare selected cluster in separate figure');
    uipushtool(th,'CData',data.ims.im_scatter,'ClickedCallBack',{@makePlotFeatures,h},'TooltipString','Plot features of selected clusters in separate figure');

    % view and execute hsplits
    uipushtool(th,'CData',data.ims.im_eye,'ClickedCallBack',{@update_comparison,h},'TooltipString','Compare selected subclusters to rest of cluster (e)','Separator','on');
    uipushtool(th,'CData',data.ims.im_go,'ClickedCallBack',{@execute_change,h,1},'Tag','executeButton','TooltipString','Execute split (x)');
    uipushtool(th,'CData',data.ims.im_go2,'ClickedCallBack',{@execute_change,h,2},'Tag','altExecuteButton','TooltipString','Merge selected and split (m)');
    
    % panel callbacks
    %figdata.heartbutton = uipushtool(th,'CData',data.ims.im_split,'ClickedCallBack',{@break_microcluster,h},'TooltipString','Split simple cluster','Enable','off');
    uipushtool(th,'CData',data.ims.im_select,'ClickedCallBack',{@select_all,h},'TooltipString','Select all panels (a)','Separator','on');
    uipushtool(th,'CData',data.ims.im_deselect,'ClickedCallBack',{@deselect_all,h},'TooltipString','Deselect all panels (d)');
    uipushtool(th,'CData',data.ims.im_hide,'ClickedCallBack',{@hide_selected,h},'TooltipString','Hide selected panels (h)');
    uipushtool(th,'CData',data.ims.im_reveal,'ClickedCallBack',{@reveal_all,h},'TooltipString','Reveal hidden panels (r)','Enable','off','Tag','reveal_button');
    uipushtool(th,'CData',data.ims.im_rearrange,'ClickedCallBack',{@rearrange_panels,h},'TooltipString','Match panel tiling to figure size (p)','Enable','on');
    
    % save spikes data
    figdata.spikes = spikes;
    figdata.original_clus = clus;
    figdata.clus   = clus;
    figdata.cutoff  = -1;
    figdata.selected = [];
    figdata.panel_ids = [];
    figdata.num_panels = 0;
    figdata.panel_list = [];
    figdata.savename = savename;
    figdata.target_fig = target_fig;
    figdata.hidden_clusters = [];

    % initialize the panels
    fw = 2;
    mw =4;
    figdata.handles.tree_panel = make_uipanel_on_grid( [fw 1 fw fw], num2str(['Click plot to select level to examine aggregation tree']), spikes.params.display, h );
    figdata.handles.cluster_panel = make_uipanel_on_grid( [fw+1 1 mw 1],'Spikes from unselected panels', spikes.params.display, h );
    figdata.handles.subcluster_panel = make_uipanel_on_grid( [fw+2 1 mw 1],'Spikes from selected panels', spikes.params.display, h );
    figdata.handles.comparison_panel = make_uipanel_on_grid( [fw fw+1 mw-fw fw],'Comparison of cluster and selection', spikes.params.display, h );
    figdata.handles.panels = [] ;
    figdata.handles.axes = [];
    figdata.handles.id = [];
            
    set(h,'UserData',figdata)

    % initialize figure panels
    display_cluster( h, clus, fw+1, figdata.handles.cluster_panel, 0 );
    init_tree(h);
    pause(.01)
    
    % set up figure
    sliderFigure( h, spikes.params.display.outer_margin, th, spikes.params.display.aspect_ratio )
    
    % set up figure
    set(h,'Name','Split Tool','NumberTitle','off')
    figure(h)
   
    
end

%%
%  callbacks for toolbar buttons
%

% save spikes object to workspace
function saveSpikes(varargin)

    h = varargin{end};
    figdata = get(h,'UserData');

    % if this tool was generated by a figure, use the savename as a callback
    if ~isempty( figdata.target_fig )
        close( h );
        pause(.01);
        figdata.savename( figdata.target_fig , figdata.spikes,figdata.original_clus );
   
    % else bring up dialog to save the data
    else
    
        varname = inputdlg('Save spike data in what variable?', 'Save to workspace', 1,{figdata.savename});
        if ~isempty(varname)
            assignin('base',varname{1},figdata.spikes)
            figdata.savename = varname{1};
            set(h,'UserData',figdata);
        end
    end
end
     
% make instance of show_clusters but just for the selected sub-clusters
function makeShowClust(varargin)

      figdata = get(varargin{end},'UserData');
      selected = figdata.selected;
      tempassigns = figdata.tempassigns;
      if isempty(selected)
                show_clusters(figdata.spikes, sort(setdiff(unique(tempassigns),0)), tempassigns);
      else
                show_clusters(figdata.spikes,selected,tempassigns);
      end
end

% make instance of compare_clusters but just for the selected sub-clusters
function makeSeparationAnalysis(varargin)

      figdata = get(varargin{end},'UserData');
      selected = figdata.selected;
            tempassigns = figdata.tempassigns;

      if isempty(selected)
          compare_clusters(figdata.spikes, sort(setdiff(unique(tempassigns),0)), tempassigns);
      else
          compare_clusters(figdata.spikes,selected, tempassigns);
      end
end

% make instance of plot_features but just for the selected sub-clusters
function makePlotFeatures(varargin)
      figdata = get(varargin{end},'UserData');
      figure
      tempassigns = figdata.tempassigns;
      selected = figdata.selected;
      if isempty(selected), selected = sort(setdiff(unique(tempassigns),0)); end
      plot_features(figdata.spikes,selected,3,0,tempassigns);
end

% updates 2 cluster comparison plots
function update_comparison( varargin)
    
    h = varargin{end};
    set(h,'Pointer','watch')
    pause(.01)
    figdata = get(h,'UserData');
    
    % get 2 populations
    selected = figdata.selected; 
    temp    = figdata.tempassigns;
  
    % throw warning if all subclusters are selected
    if length(selected) ==  length(setdiff( unique(temp), 0 ) )
                 warndlg( 'Cannot view comparison if all subclusters are selected', 'All subclusters selected');
     else
   
    
        which1  = ismember(temp, selected ); % selected clusters
        which2   = ~which1 & temp;  % unselected clusters


        % clear old plots
        delete( get(figdata.handles.comparison_panel,'Children') )
        delete( get( figdata.handles.subcluster_panel,'Children') );
        delete( get( figdata.handles.cluster_panel,'Children') );
        set( h,'UserData',figdata);

      % display the unselected panels
      if any( which2 )
         display_cluster( h,which2 , 3,figdata.handles.cluster_panel, 1 );  
      end
      
      % display the selected panels
      if any(which1)
            
            % display the selected panels
            display_cluster( h,which1 , 4,figdata.handles.subcluster_panel, 1 );

            % show the comparison, FLD an XCORR
            set( figdata.handles.comparison_panel,'Visible','off'), drawnow
             ax(1) = subplot(2,1,1,'Parent',figdata.handles.comparison_panel);
            plot_fld( figdata.spikes, which2, which1 );

             ax(2) = subplot(2,1,2,'Parent',figdata.handles.comparison_panel);
             plot_xcorr( figdata.spikes, which1, which2 );          
            set( figdata.handles.comparison_panel,'Visible','on')
             set( h,'UserData',figdata);

        end
    end
         set(h,'Pointer','arrow')

end

% perform the split operation
function execute_change(varargin)

    h = varargin{end-1};
    set(h,'Pointer','watch')
    pause(.01)
    
    mode = varargin{end};
    figdata = get(h,'UserData');
      
    allmembers = unique( figdata.spikes.info.kmeans.assigns( figdata.tempassigns ~= 0 ) );
    
    % Must select at least 1 and also unselect at least 1
    if isempty(figdata.selected)
        warndlg( 'Select one or more subclusters', 'No subclusters selected');
    elseif length(figdata.selected) ==  length(setdiff( unique(figdata.tempassigns), 0 ) )
                 warndlg( 'Cannot perform split if all subclusters are selected', 'All subclusters selected');
    else
                
        % get list of all subclusters to remove
        for j = 1:length(figdata.selected )
            subclusts = unique( figdata.spikes.info.kmeans.assigns( figdata.tempassigns == figdata.selected(j) ) );
            figdata.spikes = split_cluster(figdata.spikes,subclusts);           
        end
        
        % Mode 2 is when we are merge-and-splitting
        if mode == 2
            for j = 2:length(figdata.selected)
                figdata.spikes = merge_clusters( figdata.spikes, figdata.selected(1),figdata.selected(j) );
            end
        end
      
        % remove bad parts of map
        maps = figdata.map;
        bad = [];
        for j = 1:length(maps)
            if ismember(maps(j).to, figdata.selected), bad = [bad j]; end
        end
        figdata.map(bad) = [];
        figdata.tempassigns = make_temp_assigns(figdata.spikes,figdata.map );
   
        % if we threw away a subcluster with the original number, we need
        % to figure out what the new number is
        if ismember( figdata.clus, figdata.selected )
            available = setdiff( figdata.handles.id, figdata.selected );
            figdata.clus = available( ismember( available, figdata.spikes.assigns ) );
        end
               
        % UPDATE DISPLAY
        if ~isempty( get( figdata.handles.comparison_panel, 'Children') )   
            delete( get( figdata.handles.comparison_panel, 'Children' ) );
            delete( get( figdata.handles.subcluster_panel,'Children') );
        end
        set(h,'UserData',figdata)
        display_cluster( h, figdata.clus, 3, figdata.handles.cluster_panel, 0 );
        
        % delete selected panels
        which = find( ismember( figdata.handles.id, figdata.selected ) );
        delete( figdata.handles.panels(which) );
        figdata.handles.panels(which) = [];
        figdata.handles.axes(which) = [];
        figdata.handles.id(which) = [];
        figdata.selected = [];
        set(h,'UserData',figdata)
        update_titles(h);
        rearrange_panels(h);
        ax = draw_tree_from_map( h );
        l = line( get(ax,'XLim'), figdata.cutoff*[1 1] );
        set(l,'Color',[1 0 0 ], 'LineWidth', 2,'Tag','levelline')
 
    end
    %set(figdata.heartbutton,'Enable','off');   
    set(h,'Pointer','arrow')
    
end          

% break a minicluster
function break_microcluster( varargin )

    % get data
    h = varargin{end-1};
    which = varargin{end};
    figdata = get(h,'UserData');

    % start hour glass
    set(h,'Pointer','watch');
    pause(.01)
    
    % break the minicluster
    [figdata.spikes,newclust] = break_minicluster(figdata.spikes,which);
    figdata.clus_list = [figdata.clus_list newclust];
    
    % update tree, tempassigns, and assigns
    figdata.spikes.info.tree = [ which newclust; figdata.spikes.info.tree];
    ass = figdata.spikes.info.kmeans.assigns;
    figdata.tempassigns( ass==newclust) = newclust;
    my_clust = figdata.spikes.assigns( find( ass == which, 1) );
    figdata.spikes.assigns( ass==newclust ) = my_clust;
    
    % update target_cluster appearance
    a = find( figdata.handles.id == which );
    delete( figdata.handles.axes(a) );
    ax = axes( 'Parent',figdata.handles.panels(a) );
    plot_waveforms( figdata.spikes,figdata.tempassigns == which  );
    figdata.handles.axes(a) = ax;
    set( figdata.handles.panels(a),'ForegroundColor', [ 0 0 0] );
    
    % create a new panel for the new cluster
    if figdata.cutoff< 1
        pcolor = get( figdata.handles.panels(a),'BackgroundColor' );
        p = uipanel('Units','pixels','Visible','off','Position',get(figdata.handles.panels(1),'Position'),'Tag','waveform_panel','BackgroundColor',pcolor);
        cmenu = uicontextmenu;
        item = uimenu(cmenu, 'Label', 'Break in half', 'Callback', {@break_microcluster,h,newclust} );
        set(p,'UIContextMenu', cmenu,'Title',num2str(newclust),'ForegroundColor',[1 0 0]);             
        set( figdata.handles.panels(a),'ForegroundColor', [ 1 0 0] );
    
        a = axes('Parent',p);   
        plot_waveforms( figdata.spikes,newclust == figdata.tempassigns,3 ,1 );
        figdata.handles.panels(end+1) = p;
        figdata.handles.axes(end+1) = a;
        figdata.handles.id(end+1) = newclust;
  
        % save info
        [figdata.handles.id,a] = sort( figdata.handles.id);
        figdata.handles.axes = figdata.handles.axes(a);
        figdata.handles.panels = figdata.handles.panels(a);
        if ismember( which, figdata.selected )
            figdata.selected = [figdata.selected newclust ];
        end
    end
    
    % update map
    figdata.cutoff = figdata.cutoff + 1;
    set(h,'UserData',figdata);
        
    % update titles, panel locations, tree, break microcluster button
    rearrange_panels(h)
    ax = draw_tree_from_map( h );
    delete( findobj(ax,'Tag','levelline'));
    
    l = line( get(ax,'XLim'), figdata.cutoff*[1 1] );
    set(l,'Color',[1 0 0 ], 'LineWidth', 2,'Tag','levelline')
         
    % turn off stopwatch
    set(h,'Pointer','arrow')
    
end

% select all panels
function select_all( varargin )
    h = varargin{end};
    figdata = get(h,'UserData');
    set(figdata.handles.panels,'BackgroundColor',[1 1 1] );
    figdata.selected = setdiff( figdata.handles.id, figdata.hidden_clusters );
    set(h,'UserData',figdata)   
end

% deselect all panels
function deselect_all( varargin )
    h = varargin{end};
    figdata = get(h,'UserData');
    set(figdata.handles.panels,'BackgroundColor',figdata.spikes.params.display.label_colors(1,:) );
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
        color = figdata.spikes.params.display.label_colors(1,:);
        which = ismember( figdata.handles.id, figdata.selected );
        set(figdata.handles.panels(which),'BackgroundColor',color);
     
      figdata.selected = [];
       
       set( findobj(h,'Tag','reveal_button'), 'Enable','on');
      set(h,'UserData',figdata) 
      rearrange_panels( h );
    end
end

% reveal all hidden panels
function reveal_all(varargin)
    h = varargin{end};
    figdata = get(h,'UserData');
    figdata.hidden_clusters = [];
    set(h,'UserData',figdata); 
    rearrange_panels( h );
    set( findobj(h,'Tag','reveal_button'), 'Enable','off');

end

% rearrange panels to fit width of screen
function rearrange_panels( varargin )

    h = varargin{end};

  % make all panels and axes invisible
    figdata = get(h,'UserData');
    set( figdata.handles.panels, 'Visible','off')
    
    % get list of clusters we are supposed to show
    show_list = setdiff( figdata.handles.id, figdata.hidden_clusters );
    
    % move them to correct location
     fpos = get(h,'Position');
     d= figdata.spikes.params.display;
     r_button = findobj(h,'Tag','reset');
     if isempty(r_button), zoom_factor = 1; else, zoom_factor = get(r_button,'UserData'); end
     num_cols = max( 1, floor( (fpos(3) - d.outer_margin) / ((d.margin+d.width)/zoom_factor) ) );
    for j = 1:length(show_list)
        % where would this one go?
        row = ceil(j/num_cols);
        col = rem(j,num_cols);  if col==0,col=num_cols;end
    
       % find which one
       a = find( figdata.handles.id == show_list(j) ); 
       
       % set new position
       set(figdata.handles.panels(a), 'Position', get_pos_on_grid([4+row col 1 1], figdata.spikes.params.display, h ) );
       
    end
    % make them visible
    indices = find( ismember( figdata.handles.id, show_list ) );
    set( figdata.handles.panels(indices),'Visible','on')
    
end

% shortcut keys for accessing above callbacks
function hot_keys(varargin)

    h =varargin{end-1};
    event = varargin{end};
    switch(event.Key)
        case {'a'}, select_all(h);
        case {'d'}, deselect_all(h);
        case {'s'}, saveSpikes(h);
        case {'x'}, execute_change(h,1);
        case {'m'}, execute_change(h,2);     
        case {'e'}, update_comparison(h);
        case {'h'}, hide_selected(h);
        case {'r'}, reveal_all(h);
        case {'p'}, rearrange_panels(h);
    end
end

%%
%   tree functions
%

% initialize the aggregation tree
function init_tree(h)

   figdata = get(h,'UserData');   
   axes('Parent', figdata.handles.tree_panel)
   treeax = plot_cluster_tree( figdata.spikes, figdata.clus); 
   update_tree_level( h, treeax, figdata.cutoff );
   init_panels(h)
   set(treeax,'ButtonDownFcn',{@tree_clicked, h});
  
end

% recursively determine level of tree
function order = get_order(  tree )
   
    order = tree.order;
    if ~isempty(tree.left), order = [order get_order(tree.left) ];end
    if ~isempty(tree.right), order =[order get_order(tree.right)];end
end
        
% call back for when the tree is clicked
function tree_clicked(varargin)
 
    h = varargin{end};
     set(h,'Pointer','watch')
    pause(.01)
    
    % update tree level based on click point
    treeax = gcbo;
    level = get(treeax,'CurrentPoint');
    level = level(1,2);
    update_tree_level(h, treeax,level);
    
    % refresh all the panels
    init_panels( h );

    % save data and clean up
    figdata = get(h,'UserData');
    figdata.selected = [];
    set(h,'UserData',figdata);
    % set(figdata.heartbutton,'Enable','off');   
    set(h,'Pointer','arrow')
   
end

% update cutoff level on tree
function update_tree_level(h, ax,level)
   % this is a major change so all display information will be reset    

    axdata = get(ax,'UserData');
    figdata = get(h, 'UserData' );
    
    % if level is set to -1, use default level
    if level == -1, 
        default_to_show = figdata.spikes.params.display.initial_split_figure_panels;
        orders = sort( get_order(axdata.tree) );
        if length(orders)<default_to_show, level = 0;
        else, level = orders(1+end-default_to_show);
        end
    end
    
    % clean up all panels
    delete(figdata.handles.panels);
    figdata.handles.panels = [];
    figdata.handles.axes = [];
    figdata.handles.id = [];
    figdata.hidden_clusters = [];
    set( findobj(h,'Tag','reveal_button'), 'Enable','off');
    if ~isempty( get( figdata.handles.comparison_panel, 'Children') )
        delete( get( figdata.handles.comparison_panel, 'Children' ) );
        delete( get( figdata.handles.subcluster_panel,'Children') );
        display_cluster( h, figdata.clus, 3, figdata.handles.cluster_panel, 0 );
    end
  
    level = min( level, axdata.tree.order );
    figdata.map = get_map_from_tree( axdata.tree, level );
    figdata.tree = axdata.tree;
    figdata.cutoff = level;
    figdata.tempassigns = make_temp_assigns(figdata.spikes,figdata.map );
    set(h,'UserData',figdata)
  
    % plot the line
    set(h,'CurrentAxes',ax)
    delete( findobj(ax,'Tag','levelline'));
    l = line( get(ax,'XLim'), [level level] );
    set(l,'Color',[1 0 0 ], 'LineWidth', 2,'Tag','levelline')
    
    % blck-out irrelevant nodes
     nodes= findobj(gca,'Type','hggroup');
     set(nodes,'MarkerEdgeColor',[ 0 0 0],'MarkerFaceColor',[ 0 0 0]);
     
     % same color for display trees    
     cmap = figdata.spikes.info.kmeans.colors;
     map = figdata.map;
     for j = 1:length(map)
         c = cmap( map(j).to, :);
         set( map(j).from_dot,'MarkerEdgeColor',c,'MarkerFaceColor',c);
     end
     
    % make text
    delete( findobj(ax,'Type','Text') );
    for j = 1:length(map)
        x = get( map(j).to_dot, 'XData');
        y = get( map(j).to_dot, 'YData');
        t = text(x + .1, y, num2str(map(j).to));
        set(t,'VerticalAlignment','bottom')
    end
    
    delete( get( figdata.handles.comparison_panel,'Children' ) )
    
 
end    

% draw the tree
function ax = draw_tree_from_map( h )
    
   figdata = get(h,'UserData');
   spikes = figdata.spikes;
   which = figdata.clus;
   p = figdata.handles.tree_panel;
   map = figdata.map;
   cmap = figdata.spikes.info.kmeans.colors;

    % delete anything here
    delete( get(p,'Children') );

    % make the tree
    set(p,'Visible','off'),drawnow
    ax = axes('Parent',p);
    plot_cluster_tree( figdata.spikes, figdata.clus); 
  
    set(ax,'ButtonDownFcn',{@tree_clicked, h});
   
    axdata = get(ax,'UserData');
    map = get_map_from_tree( axdata.tree, figdata.cutoff );
  
    % black-out irrelevant nodes
     nodes= findobj(gca,'Type','hggroup');
     set(nodes,'MarkerEdgeColor',[ 0 0 0],'MarkerFaceColor',[ 0 0 0]);
     
     % same color for display trees    
     for j = 1:length(map)
         c = cmap( map(j).to, :);
         set( map(j).from_dot,'MarkerEdgeColor',c,'MarkerFaceColor',c);
     end
     
    % make text
    delete( findobj(ax,'Type','Text') );
    for j = 1:length(map)
        x = get( map(j).to_dot, 'XData');
        y = get( map(j).to_dot, 'YData');
        t = text(x + .1, y, num2str(map(j).to));
        set(t,'VerticalAlignment','bottom')
    end
    set(p,'Visible','on')

    % save map
    figdata.map = map;
    set(h,'UserData',figdata);
end
    
% make a fake assignments array from a map
function tempassigns = make_temp_assigns( spikes, map )
    
    subassigns = spikes.info.kmeans.assigns;
    assigns    = spikes.assigns;
    tempassigns = zeros( size( assigns ) );
    subcluster = unique(subassigns);
  
    for j = 1:length(map)      
        tempassigns( ismember(subassigns, map(j).from ) ) = map(j).to;
    end
end

% recursively build a map from the tree
function map = get_map_from_tree( tree, level )

    % recursively convert a graphical tree structure into a map using a cutoff level
    if tree.order <= level | isempty(tree.left)
        map.to = tree.local_id;
        map.to_dot = tree.dot;
        map.from = unique( [tree.local_id append_me(tree.left,'local_id') append_me(tree.right,'local_id') ] );
        map.from_dot = [tree.dot append_me( tree.left, 'dot' ) append_me( tree.right,'dot') ];     
    else
        map = [ get_map_from_tree( tree.left, level ), get_map_from_tree( tree.right, level ) ];
    end
end

% recursively generate depth-first ordered list of nodes of tree
function list = append_me( tree, fname )

        if isempty(tree)
            list = [];
        elseif ~isempty( tree.left )
            list = [ getfield(tree,fname) append_me(tree.left,fname) append_me(tree.right,fname) ];
        else
            list = getfield(tree,fname);
        end
end


%%
%  panel functions
%

% initialize panels
function init_panels( h )

        figdata = get(h,'UserData');
        set(h,'Pointer','watch')
        pause(.01)
    
        figdata.clus_list = setdiff( sort( unique( figdata.tempassigns ) ), 0 );
        pos =  get_pos_on_grid([5 1 1 1], figdata.spikes.params.display, h );
   
       % for each cluster
        for k = 1:length(figdata.clus_list)

            % create waveform panel
            p = uipanel('Title',num2str(figdata.clus_list(k)),'Units','pixels','Visible','off','Position',pos,'Tag','waveform_panel');
            set(p,'Title',num2str(figdata.clus_list(k)) );
            set(p,'ButtonDownFcn',{@click_panel,h, p} );
            idx = figdata.spikes.labels( figdata.spikes.labels(:,1) == figdata.clus_list(k), 2 );
            color = figdata.spikes.params.display.label_colors(1,:);
            set(p,'BackgroundColor',color);
     
            % create waveform axes
            a = axes('Parent',p);   
            plot_waveforms( figdata.spikes, figdata.clus_list(k) == figdata.tempassigns,3 ,1 );
                 
            % is this a minicluster?
            good = figdata.clus_list(k) == figdata.tempassigns;
            is_mini = length( unique( figdata.spikes.info.kmeans.assigns( good ) ) ) == 1;
            if is_mini
                set(p,'ForegroundColor',[1 0 0] ); 
                cmenu = uicontextmenu;
                item = uimenu(cmenu, 'Label', 'Break in half', 'Callback', {@break_microcluster,h,figdata.clus_list(k)} );
                set(p,'UIContextMenu', cmenu);  
            end
            
            % save info
            figdata.handles.panels(k) = p;
            figdata.handles.axes(k) = a;
            figdata.handles.id(k) = figdata.clus_list(k);

        end
      
        set( h, 'UserData',figdata);
        update_titles(h)
        rearrange_panels(h);
        set(h,'Pointer','arrow')
end

% delete selected panels
function delete_panels( ps )
    for j = 1:length(ps)
       d = get(ps(j),'UserData');
       delete( d.my_axes );
       delete(ps(j));
    end
end

% select a panel
function click_panel( varargin )
    pf = varargin{end-1};
    p  = varargin{end};
    
    pd = get(pf,'UserData');

    if ~isequal( get(pf,'SelectionType'), 'alt' )
        
        % update list
        clus = pd.handles.id( pd.handles.panels == p );
        if ismember( clus, pd.selected )
            set(p,'BackgroundColor',pd.spikes.params.display.label_colors(1,:) );
            pd.selected( pd.selected == clus ) = [];
        else        
            set(p,'BackgroundColor',[1 1 1] );
            pd.selected( end +1  ) = clus;
        end

        set(pf,'UserData',pd)
    end
end

%%
%   utility functions
%

% update cluster/subcsluter display for selected panels
function display_cluster( h, clus, row, p, issubtree )
        
    % populate cluster panel 
    if nargin < 5, issubtree = 0; end
    figdata = get(h,'UserData');
    spikes = figdata.spikes;
    set(0,'CurrentFigure',h)
    delete(get(p,'Children'))
    set(p,'Visible','off')
    drawnow
    
    % plot waveforms
    ax(1) = subplot(1,4,1,'Parent',p);
    if issubtree
        plot_waveforms(spikes,clus,3);
    else
        plot_waveforms(spikes,clus);
    end
    
    % plot detection criterion
    ax(2) = subplot(1,4,2,'Parent',p); 
    plot_detection_criterion(spikes,clus);

    % plot inter-spike interval histogram
    ax(3) = subplot(1,4,3,'Parent',p);
    plot_isi(spikes,clus);
    
    % stability plot
    ax(4) = subplot(1,4,4,'Parent',p);
    [ax(4),ax(5)] = plot_stability(spikes,clus);
   
    temp = p;
   % adjust margins -- so arbitrary! but pretty
    p = get(ax(1),'Position'); m = p(1)/3; set(ax(1),'Position',[p(1)-m p(2:4)]);
    p = get(ax(2),'Position'); set(ax(2),'Position',[p(1)-(m/3) p(2:4)]);
    p = get(ax(3),'Position'); set(ax(3),'Position',[p(1)+(m/3) p(2:4)]);
    p = get(ax(4),'Position'); set(ax([4 5]),'Position',[p(1)+m p(2:4)]);
   set(temp,'Visible','on')

end

% updates titles on each plot_waveform plot
function update_titles( h )

    % get spike train
    p = findobj(h,'Tag','waveform_panel');
    figdata = get(h,'UserData');
    ind = find(ismember(figdata.tempassigns, figdata.clus_list ));
    
    st = unwrap_time( figdata.spikes.spiketimes, figdata.spikes.trials, figdata.spikes.info.detect.dur, figdata.spikes.params.display.trial_spacing);
    
    % calcualte RPVs for entire cluster
    RPVs = sum( diff(st(ind)) < (figdata.spikes.params.refractory_period/1000) );

    % for each sub-cluster, make a title of how many RPVs are created by this sub-cluster
    for j = 1:length(p)
       a = get(p(j),'Children');
       if length(a) == 1 & isequal( get(a,'Type'), 'axes' ) % make sure "a" contains exactly 1 handle that is an axes
         clus = str2num(get(p(j),'Title'));
         ind1 = find(clus == figdata.tempassigns);
         RPVs1 = sum( diff(st(setdiff(ind,ind1))) < (figdata.spikes.params.refractory_period/1000) );
         title(a, ['Contributes ' num2str(RPVs - RPVs1) ' RPVs'] )
       end
    end

end
