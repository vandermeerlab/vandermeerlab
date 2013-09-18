function pos = get_pos_on_grid( grid_pos, params, h, impose_margins )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% get_pos_on_grid  - get 'Position' values for a new object in a sliderFigure
%
% Usage:
%    pos = get_pos_on_grid( grid_pos, params, h, impose_margins )
%
% Description:  
%    Calculate position using "pixel" coordinates of a specified object.
% This function assumes you are plotting axes or a uipanel onto a figure
% that is using the SliderFigure tool, so that it can impose a correct
% scaling and offset to the requested coordinates.  The requested 
% coordinates are given in terms of a grid position for a grid system
% that is defined by the parameters in "params".
%
% Input: 
%   grid_pos - [row column width height] of object to make
%   params   -  a display paramaters structure with fields
%          width           - number of pixels in width of plot
%          margin          - number of pixels in between plots
%          aspect_ratio    - ratio of height to width of plot
%          outer_margin    - number of pixels at figure margin
%
% Optional Inputs:
%   h              - handle to target figure (default = gcf)
%   impose_margins - whether or not to impose "margin"
%                  - typicaly, 1 for axes, 0 for panels (default =0)
%
% Output:
%  pos - [x y w h] position vector in pixels for specified graphical object
%

    % check arguments
    if nargin < 3, h = gcf; end
    if nargin < 4, impose_margins = 0; end

    
    % build as if bottom left is our top corner
    set(h,'Units','pixels')
    
    pos(1) = (grid_pos(2)-1)*(params.width+params.margin);
    pos(2) =  -grid_pos(1)* params.aspect_ratio*(params.width+params.margin);
    pos(3) =   grid_pos(3) * (params.width+params.margin) ;
    pos(4) =   grid_pos(4) * (params.width+params.margin) *params.aspect_ratio;

    % zoom it and add offsets
    z = get_zoom(h);
     pos = pos * z;

    % impose offsets
    [xoffset yoffset ] = get_offsets( params,h,z );
    pos(1:2) = pos(1:2) + [xoffset yoffset];
 
    
    % impose margins
    if impose_margins  
      pos = pos + get_zoom(h) * params.margin * .5 * [ 1 params.aspect_ratio -2 -2*params.aspect_ratio ];
    end
    
end
    
% get zoom factor which is hidden in the userdata of the reset button
function z = get_zoom(h)
   reset   = findobj(h,'Tag','reset');
   if ~isempty(reset)
       z = 1 / get( reset, 'UserData');
   else
       z = 1;
   end
end
     
% get offsets by looking at the objects in the figure
 function [x, y] = get_offsets(params, h, z)

     fig_pos = get(h,'Position');
     my_axes = findobj(h,'-depth',1,'Type','Axes');
     my_panels = findobj(h,'-depth',1,'Type','uipanel');
     xmargin = params.outer_margin;
     ymargin = params.outer_margin*params.aspect_ratio;
     
     if isempty( findobj(h,'Tag','xslider') ) | (isempty(my_axes) & isempty(my_panels) )
         x = xmargin;
         y = fig_pos(4)-ymargin;
              
     else
            % get top left position based on axes
            [axes_left,axes_top] = get_top_left( my_axes, params );
            axes_left = axes_left - .5*params.margin*z;
            axes_top  = axes_top + .5*params.margin*params.aspect_ratio*z;              
           
             % get top left position based on panels
             [panel_left,panel_top] = get_top_left( my_panels, params );
            % get overall top left
            x = min(axes_left, panel_left);
            y = max(axes_top, panel_top);
           
     end
 end 
 
% find the top left corner of a set of handles
function [l,t] = get_top_left( handles, params )
      
        if ~isempty(handles )
          pos = get(handles,'Position');
          if iscell(pos), pos = cell2mat(pos); end
          l = min(pos(:,1));
          t  = max( pos(:,2) + pos(:,4 ) );   % add left most point plus height of that object           
        else
           l = inf;  t = -inf;
        end       
end