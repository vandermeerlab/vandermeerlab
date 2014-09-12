function outlier_tool( spikes, show_clusters, h, varname,target_fig)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% outlier_tool - tool for selection of outliers from a cluster in a sorted spikes object
%
% Usage:
%      outlier_tool( spikes, show_clusters, h, varname,target_fig)
%
% Description:  
%      Tool to inspect and remove outliers from single clusters.  This tool can be called
% alone but is typically accessed through splitmerge_tool. See spike_sorting_manual.pdf 
% for instructions in how to use this tool.
%
% See ss_default_params.m as several fields in spikes.params.display affect
% the appearance of this tool.
%
%             default_figure_size          - size of figure
%             figure_font_size             - font size
%             outlier_fig_color            - color of figure background
%             label_colors                 - 1st entry sets color of panels
%             default_outlier_method       - method used to calculate distance of waveforms from center of cluster
%                                          - see plot_distances.m   
%             
% Also see sliderFigure.m since this tool makes use of sliderFigure controls.
%   
% Input: 
%   spikes         - a spikes object
%   show_clusters  - IDs for the clusters to examine
%
% Optional Input:
%   h              - handle for figure to be used for the splitmerge tool (default = gcf)
%   varname        - default name for workspace variable to save to (default = 'spikes')
%                  - this can also be a function handle to call when saving
%   target_fig     - handle for splitmerge_tool figure that called split_tool (default = [])
%
   
     % check arguments
    if ~isfield(spikes,'assigns'), error('No assignments found in spikes object.'); end
    if nargin < 2 | isequal( show_clusters,'all') 
        figdata.hidden_clusters = []; 
    else
        figdata.hidden_clusters = setdiff( unique(spikes.assigns), show_clusters );
    end
    if nargin < 3 | isempty(h), h = figure('Units','Normalized','Position',spikes.params.display.default_figure_size); end
    if nargin < 4, varname = 'spikes'; end
    figdata.mode = nargin<5;
    if ~figdata.mode,figdata.target_fig = target_fig; else, figdata.target_fig = []; end
    
    % set up figure
    clf(h,'reset');
    set(h,'defaultaxesfontsize',spikes.params.display.figure_font_size);
    set(h,'Interruptible','off','BusyAction','cancel')
    set(h,'Color',spikes.params.display.outlier_fig_color);
    set(h,'UserData',figdata,'KeyPressFcn',@hot_keys)

    % make toolbar buttons
    th = uitoolbar(h);
    data.ims = load('icons_ims');
    if figdata.mode
        tip = 'Save to workspace (s)'; icon = data.ims.im_save;
    else
       tip = 'Save to merge tool and close (s)'; icon = data.ims.im_save_and_close;
    end
    figdata.sb= uipushtool(th,'CData',icon,'ClickedCallBack',{@saveSpikes,h},'Tag','saveButton','TooltipString',tip,'Separator','on');
    if figdata.mode
      figdata.sfb = uipushtool(th,'CData',data.ims.im_savefile,'ClickedCallBack',{@saveSpikesToFile,h},'Tag','saveFileButton','TooltipString','Save to file');
      figdata.lfb = uipushtool(th,'CData',data.ims.im_loadfile,'ClickedCallBack',{@loadSpikesFromFile,h},'Tag','loadFileButton','TooltipString','Load from file'); 
    end
    figdata.gb = uipushtool(th,'CData',data.ims.im_go,'ClickedCallBack',{@execute_change,h},'Tag','executeButton','TooltipString','Remove outliers (x)','Separator','on');
     
    % save spikes data
    figdata.savename = varname;
    figdata.spikes = spikes;
    figdata.selected  = [];
    figdata.clus_list = sort(unique(spikes.assigns));
    figdata.handles.panels = [];
    figdata.method = spikes.params.display.default_outlier_method;
    figdata.cutoff = [];
    figdata.z = [];
    figdata.cutoff_line = [];
  
    % initialize the panels
    figdata.handles.scatter_panel = make_uipanel_on_grid( [2 3 2 2], 'Cluster feature plot', spikes.params.display, h );
    figdata.handles.cutoff_panel = make_uipanel_on_grid( [1 1 2 1], 'Directions:  Select cluster and then click axes to choose cutoff', spikes.params.display, h );
    figdata.handles.waveform_panel = make_uipanel_on_grid( [2 1 2 1], 'Outlier waveforms', spikes.params.display, h );
    
    % populate the panels
    set(h,'UserData',figdata)
    init_panels(h);
    figdata =get(h,'UserData');
    if ~figdata.mode
        click_panel( h,figdata.handles.panels )
    end
    
    % set up figure
    set(h,'Name','Outlier Tool','NumberTitle','off')
    sliderFigure( h, spikes.params.display.outer_margin, th, spikes.params.display.aspect_ratio )
    figure(h)
    
end

%%
%  callbacks for toolbar buttons
%

% save spikes object to workspace
function saveSpikes(varargin)
    
    h = varargin{end};
    figdata = get(h,'UserData');
    if figdata.mode % in normal mode, ask user what variable to use to save the spikes object
        varname = inputdlg('Save spike data in what variable?', 'Save to workspace', 1,{figdata.savename});
        if ~isempty(varname)
            assignin('base',varname{1},figdata.spikes)
            figdata.savename = varname;
            set(h,'UserData',figdata);
        end
    else % in the mode where outlier tool was called by splitmerge, send data back to splitmerge figure
        clus = [];
        if ~isempty(figdata.handles.id), clus = figdata.handles.id; end
      
        figdata.savename( figdata.target_fig , figdata.spikes, clus );
     
        close( h ); 
      
    end
   
end

% save spikes object to file
function saveSpikesToFile(varargin)

    h =varargin{end};    
    figdata= get(h,'UserData');
    b = figdata.sfb;
    % check if we have a previously used file name
    filename_default = get(b,'UserData');
    if isempty(filename_default)
        filename = 'spikes.mat';
    else
        filename = [filename_default.pathname filename_default.filename];
    end
    
    % ask user where to save file
    [FileName,PathName,FilterIndex] = uiputfile('*.mat','Save Spikes object',filename);
    if ~isequal(FileName,0)
        
        spikes = figdata.spikes;
        save([PathName FileName],'spikes');

        % save this name
        a.pathname = PathName;
        a.filename = FileName;
        set(b,'UserData',a);
    end
end
     
% load spikes object from file
function loadSpikesFromFile(varargin)
    
    h = varargin{end};
    
    % check if we have a default filename
    filename_default = get(findobj(h,'Tag','saveFileButton'),'UserData');
    if isempty(filename_default)
        filename = [pwd '\spikes.mat'];
    else
        filename = [filename_default.pathname filename_default.filename];
    end
    
    % ask user where file is
    [filename,pathname] = uigetfile('*.mat','Load spikes object from file.',filename);
   
    if ~isequal(filename,0)
      
      % load data
      a = load([pathname filename]);  
      names = fieldnames(a);
      figdata = get(h,'UserData');
      spikes = getfield(a,names{1});
      
      % do a hard reset
      outlier_tool( spikes, unique(spikes.assigns), h, figdata.savename)
    
    end
end
 
% remove outliers over set threshold
function execute_change(varargin)

    % start stopwatch
    h = varargin{end};
    figdata = get(h,'UserData');
    set(h,'Pointer','watch'),pause(.01)
 
    % a cluster must be selected
    if isempty(figdata.selected)
        warndlg( 'Select a cluster to remove outliers from.', 'No cluster selected');
    else
    
        % call remove outliers
        indices = find( figdata.spikes.assigns == figdata.selected );
        baddies = find( figdata.z > figdata.cutoff );
        my_list = zeros( [1 length( figdata.spikes.assigns ) ] );
        my_list(indices(baddies) ) = 1;
        my_subclusters = unique( figdata.spikes.info.kmeans.assigns( indices ) );
        figdata.spikes = remove_outliers(figdata.spikes, my_list  );
        
        % check whether the user completely emptied out the cluster (bad user)
        a = find( ismember( figdata.spikes.assigns, my_subclusters ), 1);
        old_selected = figdata.selected;
           
        figdata.selected = figdata.spikes.assigns(a);
        if isempty(figdata.selected), figdata.selected = NaN; end
        
       set(h,'UserData',figdata)
       if sum( figdata.spikes.assigns == figdata.selected ) == 0
        
            % clear everything
            which = find(figdata.selected == figdata.handles.id);
            figdata.selected = [];
            delete(figdata.handles.panels(which));
            figdata.handles.id(which) = [];
            figdata.handles.panels(which) = [];
            figdata.handles.axes(which) = [];
            
            % delete the other panels
            delete( get( figdata.handles.scatter_panel,'Children') );
            delete( get( figdata.handles.cutoff_panel,'Children') );
            delete( get( figdata.handles.waveform_panel,'Children') );

            rearrange_panels(h)
            set(h,'UserData',figdata)
            
       else % update display for fewer events
         
            % find axes to update
            which =  find( figdata.handles.id == old_selected );
            ax = figdata.handles.axes(which);
            set(h,'UserData',figdata);
            set(figdata.handles.panels(which),'Title',num2str(figdata.selected) );
            
            % update my panel, zvalue_histogram, outlier waveforms, scatter
            ud = get(ax,'UserData');      
            colormode = ud.colormode;
            set(h,'CurrentAxes',ax)
            plot_waveforms( figdata.spikes, figdata.selected, colormode,0);
            make_zvalue_histogram(h)
            update_waveforms( h)
            update_scatter( h )
            
        end
        
    end
    set(h,'Pointer','arrow')
 
end      

% shortcut keys for accessing these callbacks
function hot_keys(varargin)

    h =varargin{end-1};
    event = varargin{end};
    switch(event.Key)
        case {'s'}, saveSpikes(h);
        case {'x'}, execute_change(h);
    end
end
   
%% 
%   panel functions
%

% initialize all panels
function init_panels( h )

    % start stopwatch
    figdata = get(h,'UserData');
    set(h,'Pointer','watch')
    pause(.01)
    
    % for each cluster
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
        set(p,'UIContextMenu',get_panel_menu( figdata.spikes.params.display.label_categories, idx, p,h ) );  
        
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

% move them to proper location
function rearrange_panels( h )

    % make all panels and axes invisible
    figdata = get(h,'UserData');
    set( figdata.handles.panels, 'Visible','off')
    show_list = setdiff(figdata.spikes.assigns, figdata.hidden_clusters);
    
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
       set(figdata.handles.panels(a), 'Position', get_pos_on_grid([2+row col 1 1], figdata.spikes.params.display, h ) );
       
    end
    
    % make them visible
    indices = find( ismember( figdata.handles.id, show_list ) );
    set( figdata.handles.panels(indices),'Visible','on')
    
end

% select/deselect a panel
function click_panel( varargin )
    h = varargin{end-1};
    
    % must be a left-click
    if isequal( get(h,'SelectionType'), 'normal' )
    
       % start stopwatch
       set(h,'Pointer','watch'),pause(.01)
        p  = varargin{end};

        figdata = get(h,'UserData');
        clus = str2num(get(p,'Title'));
        
        % check for (impossible?) case of no spikes
        if sum( figdata.spikes.assigns == clus ) == 0
              warndlg( 'Cannot select cluster.  All spikes were previously labeled outliers.','Empty cluster.');
        else

            % selection is disabled in splitermerge_tool mode
            if figdata.mode
              delete( get( figdata.handles.scatter_panel,'Children') );
              delete( get( figdata.handles.cutoff_panel,'Children') );
              delete( get( figdata.handles.waveform_panel,'Children') );
            end

            % unselecting a cluster
            if figdata.selected == clus & figdata.mode
                figdata.selected = [];
                which = find( figdata.spikes.labels(:,1) == clus );
                cat   = figdata.spikes.labels(which,2);
                set(p,'BackgroundColor', figdata.spikes.params.display.label_colors( cat, : ) );
                set(h,'UserData',figdata)

            % selecting a cluster
            elseif ~isequal( figdata.selected, clus)

                %deselect the old choice
                if ~isempty( figdata.selected )             
                    which = find( figdata.spikes.labels(:,1) == figdata.selected );
                     cat   = figdata.spikes.labels(which,2);
                    set(figdata.selectedp,'BackgroundColor',figdata.spikes.params.display.label_colors( cat, : ) ); 
                end
                
                % highlight the new choice
                figdata.selected = clus;
                figdata.selectedp = p;
                set(p,'BackgroundColor', [1 1 1] );
                figdata.show_axes(1) =   axes('Parent',figdata.handles.cutoff_panel );
                figdata.show_axes(2) =   axes('Parent',figdata.handles.waveform_panel);
                figdata.show_axes(3) =   axes('Parent',figdata.handles.scatter_panel);

                % update all panels
                set(h,'UserData',figdata);
                make_zvalue_histogram( h )
                update_waveforms( h)
                update_scatter( h )
            end
         end
           set(h,'Pointer','arrow')
        end
end

% delete panels
function delete_panels( ps )
    for j = 1:length(ps)
       d = get(ps(j),'UserData');
       delete( d.my_axes );
       delete(ps(j));
    end
end

% change label for a cluster
function change_label(varargin)
    
    % user can change cluster labels via a context menu
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

% context menu for a panel
function cmenu = get_panel_menu( categories, me,p,h)
% context menu for changing cluster label
    cmenu = uicontextmenu;
    clus = str2num( get(p,'Title' ) );

    for j = 1:length(categories)
        d(j) = uimenu(cmenu, 'Label', categories{j}, 'Callback', {@change_label, p, j,h});
        if j == me, 
            set( d(j), 'Checked','on','Enable','off'); 
        end
    end
end

%%
%  callback for distances histogram
%

% populate histogram of distances
function make_zvalue_histogram( h )

    % get figure info
    figdata = get(h,'UserData');
    
    % plot the histogram
    ax = figdata.show_axes(1);
    set(h,'CurrentAxes',ax)
    [figdata.z,dof] = plot_distances( figdata.spikes, figdata.selected, figdata.method );
    
    % add a context menu to change distance calculation mode
    cmenu = uicontextmenu;
    uimenu(cmenu, 'Label', 'Use cluster stats', 'Callback', {@change_mode, h, 1},'Checked',bool2word(figdata.method ==1) );
    uimenu(cmenu, 'Label', 'Use noise stats', 'Callback', {@change_mode, h, 2},'Checked',bool2word(figdata.method ==2) );
    set(figdata.show_axes,'UIContextMenu', cmenu )
 
    % place the vertical line representing the threshold for whether a spike is a cutoff
    temp = .5/length(figdata.z); % make odds of there being any data this extreme less than 50%
    figdata.cutoff =  min( chi2inv(1-temp, dof ), max( figdata.z) );
    l = line(figdata.cutoff*[1 1],get(ax,'YLim'));
    set(l,'Color',[1 0 0],'LineWidth',2,'LineStyle','--','Tag','cutoff')
    figdata.cutoff_line = l;
    set(ax,'ButtonDownFcn', {@update_cutoff, h} );
    set(h,'UserData',figdata);
 
    % get indices of RPVs and scatter plot them
    select = get_spike_indices(figdata.spikes, figdata.selected );
    spiketimes =  sort( figdata.spikes.unwrapped_times(select) );
    which = find( diff(spiketimes) < figdata.spikes.params.refractory_period/1000 );
    which = unique( [which which+1 ] );
    hold on
    s = scatter( figdata.z(which), mean(get(gca,'YLim'))*ones(size(figdata.z(which))),'k' );
    set(s,'Marker','x')
    hold off
end

% callback for when new threshold is selected
function update_cutoff( varargin )

    h = varargin{end};
    % make sure this is a left-click
    if isequal( get(h,'SelectionType'), 'normal' )

        % start stopwatch
        set(h,'Pointer','watch')
         pause(.01)

         %update cutoff
        figdata = get(h,'UserData');
        ax = figdata.show_axes(1);
        cp = get(ax,'CurrentPoint');
        figdata.cutoff = cp(1);
        set(figdata.cutoff_line,'XData',figdata.cutoff*[1 1]);
        set(h,'UserData',figdata);
        
        % update other plots
        update_waveforms( h)
        update_scatter(h)
        set(h,'Pointer','arrow')
    end   
 end

% change to covariance matrix used by histogram, accessed by context menu
function change_mode( varargin )

    % get new method
    h = varargin{ end - 1};
    figdata = get(h,'UserData');
    method = varargin{end};
    
    % only act if method was changed
    if ~isequal( method, figdata.method )
 
       % start stopwatch
       set(h,'Pointer','watch')
       pause(.01)
 
       % save method
       figdata.method =  method;
       set(h,'UserData',figdata);
    
       % update plots
        make_zvalue_histogram( h )
        update_waveforms( h)
        update_scatter(h)
        set(h,'Pointer','arrow')
        
    end
   
end

%%
%   utility functions
%            

% updates waveform panel
function update_waveforms( h)

    % get indices of waveforms and outliers
    figdata = get(h,'UserData');
    ax = figdata.show_axes(2);
    set(h,'CurrentAxes',ax);
    cla;
    indices = find( figdata.spikes.assigns == figdata.selected );
    baddies = find( figdata.z > figdata.cutoff );
       
    % get the waveforms
    w = figdata.spikes.waveforms( indices(baddies), : );
    figdata.spikes.params.display.default_waveformmode = 1;
    which = ismember( 1:length(figdata.spikes.assigns), indices(baddies));
    if length(baddies) > 0
     plot_waveforms(figdata.spikes,which,1); 
     delete( get(ax,'UIContextMenu'))
    end
    
    % plot mean waveform
    which = find( figdata.spikes.assigns == figdata.selected );
    num_samples = size( figdata.spikes.waveforms(:,:),2);
    l = line( 1:num_samples, mean( figdata.spikes.waveforms( which ,:)));
    set(l,'Color',[ 0 0 0],'LineWidth',2.5)
   
    % set labels
    xlabel( 'Sample');
    ylabel('Value');
    set(ax,'XLim',[1 num_samples])  
    title(['Outliers for cluster #' num2str(figdata.selected) ' (' num2str(length(baddies)) ' out of ' num2str( length(figdata.z) ) ')'])

end

% updates scatter plot
function update_scatter( h )

    % get figure and axes data
    figdata = get(h,'UserData');
    ax = figdata.show_axes(3);
    set(h,'CurrentAxes',ax);
    axdata = get(ax,'UserData');
    
    % get regular waveforms and outliers
    indices = find( figdata.spikes.assigns == figdata.selected );
    baddies = find( figdata.z > figdata.cutoff );
    s = figdata.spikes;
    c = max(s.info.kmeans.assigns) + 1;
    s.assigns( s.assigns ~= figdata.selected ) = 0;   
    s.assigns( indices(baddies) ) = c;
    s.info.kmeans.assigns( indices(baddies) ) = c;
    s.info.kmeans.colors(c,:) = 0;
   
    % update defaults if present already
     if ~isempty(axdata)
          s.params.display.xchoice = axdata.xchoice;
          s.params.display.xparam = axdata.xparam;
          s.params.display.ychoice = axdata.ychoice;
          s.params.display.yparam = axdata.yparam; 
    end
   
    % plot the cluster with outliers in black
    plot_features(s, [c figdata.selected], 3, 0 );
    cmenu = get(ax,'UIContextMenu');
    kids  = get( cmenu, 'Children');
    delete( setdiff( kids, findobj(cmenu,'Tag','showdensityoption')))
       
end

% converts 0/1 to 'off'/'on'
function word = bool2word( b )
    if b, word = 'on'; else, word = 'off'; end
end
