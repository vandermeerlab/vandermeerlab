function panel = make_uipanel_on_grid( panel_pos, my_title, params, h,visibility)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% make_uipanel_on_grid - returns uipanel object with specified position
%
% Usage:
%   panel = make_uipanel_on_grid( panel_pos, my_title, params, h,visibility)
%
% Description:  
%   This is a utility function that generates a blank uipanel object at 
% a specified position.  The position is defined by both the sliderFigure,
% which can scale panels placed in the figure, and by the display params 
% passed in as input.  The logic of panel placement is in the function 
% get_pos_on_grid. 
%
% Input: 
%   axes_pos   - [row column width height] location of uipanel to create
%              - see get_pos_on_grid.m for how to interpret these numbers
%   my_title   - title of the uipanel
%   params     - params object
%              - see get_pos_on_grid.m for specification of parameters
%
% Optional inputs:
%   h          - handle to target figure (default is gcf)
%   visibility - 'on' or 'off', should the panel be visible? (default is 'on')
%
% Output:
%   panel      - handle to panel at the specified location
%

    % by default, the target figure is the current figure and visibility is
    % turned on
    if nargin < 4, h = gcf; end
    if nargin < 5, visibility = 'on'; end
    
    % get position on grid
    pos = get_pos_on_grid( panel_pos, params, h );

    % make panel
    set(0,'CurrentFigure',h)
    panel = uipanel('Title',my_title,'Units','pixels','Position',pos,'Visible',visibility);
    
    % our panels will keep track of their axes
    data.my_axes = [];
    set(panel,'UserData',data)
    
    % reset the slider tool
    reset_slider(h);

    
    