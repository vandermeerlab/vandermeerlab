function sliderFigure( h, margin, th, yar )
% SliderFigure tool by Hill DN - 07/12/2010
%
% sliderFigure - Allows user to scale and scroll through all axes on a figure
%
% Usage:
%      sliderFigure( h, margin, toolbar, yar )
%
% Description:  
%   Call sliderFigure to add scroll bars and a toolbar with zoom control 
%  buttons.  If additional plots are added to the figure later, call 
%  reset_slider.m to refreseh the state of the sliderFigure.
%
%   The function of the zoom control buttons from left to right are:
%   (1) zoom in (10%)
%   (2) zoom out (10%)
%   (3) fit all objects to the width of the current figure
%   (4) fit all objects to the height of the current figure
%   (5) fit all objects to both the width and height of the current figure
%   (6) reset zoom to 100% 
%  
% Notes:
%   (1) sliderFigure only supports figures composed of uipanels and axes.
%   (2) All 'units' fields for uipanels and axes are changed to 'pixels'
%   (3) All axes/uipanels that have coordinates in a visible location will
%       have their "visibility" turned "on".  All others will be set to "off".
%   (4) To have the slider update the display in real time, it needs to be
%       seen as the current object.  Otherwise the display will not update
%       until you release the mouse button.  You may need to click the 
%       slider twice to make it the current object.
%
% Optional Inputs: 
%     h            - handle to target figure (default is to use gcf)
%     margin       - width of border to impose in pixels (default = 60)
%     toolbar      - handle to uitoolbar object to append sliderFigure buttons to (default is a new uitoolbar)
%     aspect_ratio - ratio of width of vertical border to width of horizontal border (default = 3/4)
%

% Check arguments
if nargin < 1, h =gcf;       end
if nargin < 2, margin = 60;  end
if nargin < 4, yar = 3/4;   end

% move focus to current figure
set(0,'CurrentFigure',h);
if isempty(findobj(h,'Tag','xslider'))

    % create controls
    xslider = uicontrol('style','slider','units','normalized', ...
       'value', 0, ...
       'callback',{@slider_callback,h}, ...
       'userdata',margin, ...
        'Tag','xslider');

    yslider = uicontrol('style','slider','units','normalized', ...
       'value', 1, ...
       'callback',{@slider_callback,h}, ...
       'userdata', yar*margin,...
             'Tag', 'yslider' );

    corner_button = uicontrol('Style', 'pushbutton', 'String', '', 'Tag','corner','Enable','off');
 
    % create tool bar with its various icons
    if nargin < 3,th = uitoolbar(h); sep_val ='off'; else, sep_val='on'; end
    
    ims = load('slider_ims'); 
    uipushtool(th,'CData',ims.im_plus,'ClickedCallBack',{@zoom_it,h},'userdata',[1.25],'TooltipString','Zoom in','Separator',sep_val);
    uipushtool(th,'CData',ims.im_minus,'ClickedCallBack',{@zoom_it,h},'userdata',[.8 ],'TooltipString','Zoom out') ; 
    uipushtool(th,'CData',ims.im_expand,'ClickedCallBack',{@fit_callback,[1 1],h},'userdata',1,'Tag','fittoscreen','TooltipString','Fit to screen');
    uipushtool(th,'CData',ims.im_height,'ClickedCallBack',{@fit_callback,[0 1],h},'userdata',1,'Tag','fitheight','TooltipString','Fit height');
    uipushtool(th,'CData',ims.im_width,'ClickedCallBack',{@fit_callback,[1 0],h},'userdata',1,'Tag','fitwidth','TooltipString','Fit width to screen');
    uipushtool(th,'CData',ims.im_reset,'ClickedCallBack',{@reset_callback,h},'userdata',1,'Tag','reset','TooltipString','Original size');

    % set callbacks
    set(h,'ResizeFcn',{@resize_callback,h})
    set(h,'WindowButtonMotionFcn', {@sliding_callback,h});
    
end

% if the figure is not empty
if ~isempty( get_scalable_objects( h) ) 
 % call the callbacks to initialize appearance
 resize_callback(h)
end

% set current object to the corner button
 set(h,'CurrentObject',findobj(h,'Tag','corner'),'Toolbar','figure')
 
% call back for zoom in/out buttons.
% the zoom factor is stored in the calling button's userData
function zoom_it(varargin)
    
    h = varargin{end};
    
    % find the middle pixel of the figure
    set(h,'Units','Normalized'); fig_pos = get(h,'Position');
    middle = fig_pos(3:4)/2;
    
    % get other user variables
    button = get(gcbo);
    factor = button.UserData;
    my_objects = get_scalable_objects( h);
    
    % save total zoom factor in the reset button
    reset_button = findobj( h,'Tag','reset');
    x = get(reset_button,'userdata');
    set(reset_button,'userdata',x/factor);
 
    % reset axes
    for j = 1:length(my_objects)
        old_pos = get(my_objects(j),'Position');
        new_pos = [  middle + factor*(old_pos(1:2)-middle) factor*old_pos(3:4) ]; 
        set(my_objects(j),'Position', new_pos )
    end
    
    % reset sliders
    resize_callback(h)
    slider_callback(h)
    
% Callback for the fit-zoom buttons 
% Zoom factor needs to be calculated
function fit_callback(varargin)

% get all the objects of interest
h = varargin{end};
options = varargin{end-1};

xslider = findobj(h,'Tag','xslider');
yslider = findobj(h,'Tag','yslider');

% reset sliders
set(xslider,'Value',0,'Visible','off');
set(yslider,'Value',1,'Visible','off');

% what is the size of my axes
 xmargin = get(xslider,'UserData');
 ymargin = get(yslider,'UserData');
 my_objects = get_scalable_objects( h);
 pos  = get(my_objects,'Position') ; if iscell(pos), pos = cell2mat(pos); end        
 height   =  max( pos(:,2) + pos(:,4) ) - min( pos(:,2) ); 
 width   =  max( pos(:,1) + pos(:,3) ) - min( pos(:,1) );
 fig_pos     = get(h,'Position');
   
% what is the needed zoom factor
zooms = (fig_pos(3:4) - 2*[xmargin ymargin] - [20 20]) ./ [width height];

if ~options(1), zoom_factor = zooms(2);
elseif ~options(2), zoom_factor = zooms(1);
else, zoom_factor = min(zooms); end

% set that zoom
set( gcbo,'UserData', max( zoom_factor, .01)*.99);

zoom_it(h)
  
% callback for reset button. 
function reset_callback(varargin)

    % get all the objects of interest
    h = varargin{end};
    xslider = findobj(h,'Tag','xslider');
    yslider = findobj(h,'Tag','yslider');

    % reset sliders
    set(xslider,'Value',0);
    set(yslider,'Value',1);
    zoom_it(h)
         
% callback for resizing axes and uipanels 
function resize_callback(varargin)

    h = varargin{end};
    barwidth = 20; % pixels
    
    % get all the objects of interest
    xslider = findobj(h,'Tag','xslider');
    yslider = findobj(h,'Tag','yslider');
    corner_button = findobj( h,'Tag','corner');
    my_objects = get_scalable_objects( h);
   
    if ~isempty( my_objects )
    
        set( [my_objects h xslider yslider corner_button], 'Units','pixels')

        % get some numbers
        xmargin  = get(xslider,'userdata');
        ymargin  = get(yslider,'userdata');        
        pos  = get(my_objects,'Position') ; if iscell(pos), pos = cell2mat(pos); end        
        height   =  max( pos(:,2) + pos(:,4) ) - min( pos(:,2) ) + 2*ymargin; 
             
        width   =  max( pos(:,1) + pos(:,3) ) - min( pos(:,1) ) + 2*xmargin + 20;
        fig_pos     = get(h,'Position');

        % set visibility
        if fig_pos(4) > height, set(yslider,'Visible','off','Value',1); else, set(yslider,'Visible','on'); end
        if fig_pos(3) > width, set(xslider,'Visible','off','Value',0); else, set(xslider,'Visible','on'); end
        if fig_pos(4) > height & fig_pos(3) > width, set(corner_button,'Visible','off'); else, set(corner_button,'Visible','on'); end

        % set bar position
        set([xslider yslider],'Units','pixels')
        set(xslider,'Position',[ 0 0 fig_pos(3)-barwidth barwidth] );
        set(yslider,'Position',[  fig_pos(3)+1-barwidth barwidth barwidth fig_pos(4)-barwidth] );
        set(corner_button,'Position', [ fig_pos(3)+1-barwidth 0 barwidth barwidth]);
       
        % set slider width
        set(xslider, 'Units','Normalized','SliderStep', [0.01 fig_pos(3) / width ] );
        set(yslider, 'Units','Normalized','SliderStep', [0.01 fig_pos(4) / height ] );
        set([xslider yslider],'Units','pixels');
    end
        
    % finally, need to call slider callback to add offset to zoom
    slider_callback(h);

% callback to add offset position to figure view
function slider_callback(varargin)
     
    h = varargin{end};
    % get all the objects of interest
    xslider = findobj(h,'Tag','xslider');
    yslider = findobj(h,'Tag','yslider');
    my_objects = get_scalable_objects( h);
  
    % if there are any axes or uipanels
    if ~isempty(my_objects )
    
        set( [my_objects h xslider yslider], 'Units','pixels')
        
        % get some numbers
        xmargin  = get(xslider,'userdata');
        ymargin  = get(yslider,'userdata');
        
        pos  = get(my_objects,'Position') ; if iscell(pos), pos = cell2mat(pos); end        
        height   =  max( pos(:,2) + pos(:,4) ) - min( pos(:,2) ) + 2*ymargin + 20; 
        width   =  max( pos(:,1) + pos(:,3) ) - min( pos(:,1) ) + 2*xmargin + 20;
        fig_pos     = get(h,'Position');

        % get current location
        current_top  = max( pos(:,2) + pos(:,4) );
        current_left    = min( pos(:,1) );

        % figure out the  target position of everything
        xval = get(xslider,'Value');  yval = get(yslider,'Value');
        new_top = round(height - ymargin + yval*(fig_pos(4) - height) ); 
        new_left   = round( xmargin + min( xval*(fig_pos(3) - width), 0 ));
           
        % update axes
        for j = 1:length(my_objects)
              x = get(my_objects(j),'Position');
              new_pos = [x(1)+new_left-current_left  x(2)+new_top-current_top  x(3) x(4)];
              visibility = visibilty_from_pos(new_pos,fig_pos); % turn off visibility for objects not on the screen
              set(my_objects(j),'Visible',visibility);
              set(my_objects(j),'Position', new_pos)
        end
    end
       
        set(h, 'CurrentObject', xslider)
 
%if an object is not on the screen at all, then set its visibility to off
function visibility = visibilty_from_pos(new_pos,fig_pos)

    cond1 = new_pos(1) > fig_pos(3);
    cond2 = new_pos(1) + new_pos(3) < 1;
    cond3 = new_pos(2) > fig_pos(4);
    cond4 = new_pos(2) + new_pos(4) < 1;
    if cond1 | cond2 | cond3 | cond4 
        visibility = 'off';
    else
        visibility = 'on';
    end

% if current object is a slider, only then call the slider callback
function sliding_callback(varargin)
  
   if isequal( get(gco,'Tag'), 'xslider') | isequal( get(gco,'Tag'), 'yslider') 
           slider_callback(varargin{end})       
   end
   
   % get any axea and uipanels in the plot
function objects = get_scalable_objects( h)
  
    objects =  [ findobj( get(h,'Children'), 'flat','type','axes','-and')' findobj(get(h,'Children'),'flat','type','uipanel','-and')' ];
    