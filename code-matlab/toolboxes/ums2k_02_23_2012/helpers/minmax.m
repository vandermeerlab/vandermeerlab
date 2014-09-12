function [extrema,inds] = minmax(X)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/08/2010
%
% minmax - Simultaneous overall smallest and largest components. 
%
% Usage:
%    [extrema,inds] = minmax(X)
%
% Description:
%   Y = MINMAX(X) returns the minimum and maximum values, of the array X
%   such that Y(1) = MIN(X) and Y(2) = MAX(X).  For N-D arrays X,
%   MINMAX(X) is equivalent to MINMAX(X(:)).
%
%   [Y,I] = MINMAX(X) also returns the linear indices of the extrema such
%   that, Y(1) == X(I(1)) and Y(2) == X(I(2)).  When X has more than one
%   extremal element, the index of the first is returned.
%
% Input: 
%   X  - 1 dimensional array.  NaN's are ignored.
%
% Output:
%   extrema - [min(X) max(x)]
%   inds    - location of extrema within X
%

[extrema(1),inds(1)] = min(X(:));
[extrema(2),inds(2)] = max(X(:));

