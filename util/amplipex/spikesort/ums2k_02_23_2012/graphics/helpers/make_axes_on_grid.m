function ax = make_axes_on_grid( axes_pos, params, h, visibility)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% make_axes_on_grid - returns axes object with specified position
%
% Usage:
%   ax = make_axes_on_grid( axes_pos, params, h, visibility)
%
% Description:  
%   This is a utility function that generates a blank axes object at 
% a specified position.  The position is defined by both the sliderFigure,
% which can scale axes placed in the figure, and by the display params 
% passed in as input.  The logic of axes placement is in the function 
% get_pos_on_grid. 
%
% Input: 
%   axes_pos   - [row column width height] location of axes to create
%              - see get_pos_on_grid.m for how to interpret these numbers
%   params     - params object
%              - see get_pos_on_grid.m for specification of parameters
%
% Optional inputs:
%   h          - handle to target figure (default is gcf)
%   visibility - 'on' or 'off', should the axes be visible? (default is 'on')
%
% Output:
%   ax         - handle to axes at the specified location
%

    % by default, the target figure is the current figure and visibility is
    % turned on
    if nargin < 3, h = gcf; end
    if nargin < 4, visibility = 'on'; end
    
    % get grid position
    pos = get_pos_on_grid( axes_pos, params, h,1 );
        
    % build axes using the current figure and place the axes on top
    set(0,'CurrentFigure',h);
    ax = axes('Units','pixels','Position',pos,'Visible',visibility);
    uistack(ax, 'top');
    set(gcf,'CurrentAxes',ax)
    
    % reset the slider tool
    reset_slider(h);
    