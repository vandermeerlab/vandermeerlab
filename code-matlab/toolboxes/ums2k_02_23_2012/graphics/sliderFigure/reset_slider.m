function reset_slider(h)
% sliderFigure tool by Hill DN  - 07/12/2010
%
% reset_slider - reset controls on a sliderFigure figure
%
% Usage:
%      reset_slider(h)
%
% Description:  
%   Utility to reset the controls on a sliderFigure.  This should be called
% if sliderFigure gets into a bad state, or whenever new axes are drawn to
% a figure after sliderFigure is called.
%
% Optional input: 
%   h  - figure handle, must already have the sliderFigure tool running
%      - defaults to current figure
%

% check argument
if nargin<1, h =gcf; end

% grab the sliders
xslider = findobj( h,'Tag','xslider');
yslider = findobj( h,'Tag','yslider');

% do nothing if they don't exist
if ~isempty(xslider)

    %
    % get handles and variables
    %
    axes_pos = get( findobj(h,'Type','axes'),'Position' );
    if iscell(axes_pos), axes_pos = cell2mat(axes_pos); end
    panel_pos = get( findobj(h,'Type','uipanel'),'Position');
    if iscell(panel_pos), panel_pos = cell2mat(panel_pos); end
    pos = [axes_pos; panel_pos];
    fig_pos = get(h,'Position');
    xmargin = get(xslider,'UserData');
    ymargin = get(yslider,'UserData');
    
    
    % find true width and height of content of figure
    height = 2*ymargin + max(pos(:,2) + pos(:,4)) + -min(pos(:,2)) + 20;
    width  = 2*xmargin + max(pos(:,1)+pos(:,3)) - min(pos(:,1)) + 20;
 
    % find current position
    curx  = xmargin - min(pos(:,1) );
    cury  = max(pos(:,2)+pos(:,4)) + ymargin - fig_pos(4) ;

    % update slider position
    set(xslider, 'Value',  max(0, curx / ( width - fig_pos(3) ) ) );
    set(yslider, 'Value',  min( max(0, 1-(cury / ( height - fig_pos(4) ) )), 1));

    % update width of the slider itself
    a = get(xslider,'SliderStep');
    set(xslider,'SliderStep',[a(1) max( fig_pos(3)/( width - fig_pos(3)), .01 )])
    set(yslider,'SliderStep',[a(1) max( fig_pos(4)/( height - fig_pos(4)),.01 )])
  
end

