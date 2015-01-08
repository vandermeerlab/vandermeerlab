function nan_imagesc_ec(data, varargin)    
%  nan_imagesc_ec: this is a walkaround to get nans to appear as white in
%  2D grid as if it was imagesc.  This is an issue with imagesc due to the
%  two types of renderers ('painters' and 'opengl') either failing to
%  produce the nans properly or causing problems with the text on top of a
%  figure.  
%
% this will simply create a 3d surfplot and then view it from above. 
% Credit belongs to "Richard" on matlabcentral http://www.mathworks.com/matlabcentral/answers/9394#answer_13000
%
% modified and converted into a fcn by EC 2015-01-07
%%
% set defaults, 
nan_colour = 'k'; 
extract_varargin
%%
    A = data;
    % Here the NaNs show up as the colour specified in "nan_colour"
    surf(0.5:(size(A, 2)+0.5), 0.5:(size(A, 1)+0.5), ...
        zeros(size(A)+1), A, 'edgecolor', 'none');
    set(gca, 'color', nan_colour) 
    axis('tight');
    view(2);
    grid('off');
    box('on');
    set(gca, 'YDir', 'reverse');