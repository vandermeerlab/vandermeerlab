function [counts,x_inds,y_inds] = histxy(x, y, D, filt)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% histxy - 2D histogram of joint x and y values
%
% Usage:
%    [counts,x_inds,y_inds] = histxy(x, y, D, filt)
%
% Description:  
%  Returns a matrix containing counts for how many data points are found
% within a 2-dimensional bin.  The number of bins can be specified by 
% the user.  The matrix can also be filtered using a local average.
%
% Input: 
%   x  - array of scalar "x-values"
%   y  - array of corresponding "y-values"
%
% Optional inputs:
%    D    - [n_x n_y] number of bins in each dimension (default is [25 25])
%         - if only 1 value is specified, both dimensions will use that number of bins
%    filt - size of square uniform averaging matrix (default is no filtering)
%
% Output:
%   counts    - matrix of size D containing counts for each bin
%   x_inds    - location of x_vals bins
%   y_inds    - location of y_vals bins
%
  
    if nargin < 3, D = [25 25];end
    if nargin < 4, filt = []; end
    if length(D) == 1, D = [D D]; end

    % convert values to (x,y) coordinates and then indices
    [x,oldminx,oldmaxx] = rescale(x,1,D(1));  x = round(x);
    [y,oldminy,oldmaxy] = rescale(y,1,D(2));  y = round(y);
    inds = x + (y-1)* D(1);
 
    % fill in the counts based on indices
    counts = zeros( D );
    for j = 1:length(inds)
        counts(inds(j)) = counts(inds(j))+1;
    end
    
    % save output of bin coordinates
    x_inds = [oldminx oldmaxx];
    y_inds = [oldminy oldmaxy];
    
    % smooth 2D-histogram
    if ~isempty(filt)
       counts = conv2( counts, ones(filt), 'same'); 
    end
    
end
        
