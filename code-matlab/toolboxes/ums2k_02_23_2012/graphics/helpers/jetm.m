function map = jetm(n)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% jetm -  Muted jet colormap.
%
% Usage:
%      map = jetm(n)
%
% Description: 
%   JETM(N), a variant of JET(N), is an N-by-3 matrix containing a muted
%   version of the JET colormap.  JETM with no arguments is the same
%   length as the current colormap.
%

if (nargin < 1), n = size(colormap,1);  end;

% This colormap was constructed manually by altering the JET colormap,
%  which is visually a little too bright and suffers from light blue and
%  yellow bands.  To obtain a muted version,
%    (1) map a jet colormap to HSV color space
%    (2) linearly rescale the saturation values from the
%               range [0.5 1.0] => [0.25 0.65] (softens the colors)
%    (3) linearly rescale the color value numbers from the
%               range [0.5 1.0] => [0.8 1.0]  (keeps the colors bright)
%    (4) use the colormapeditor to better space out the light blue/yellow
%               color points
%    (6) remap from HSV to RGB color space
%    (5) symmetrize the map by enforcing the red/blue values to be
%               inverted copies of one another and the green values
%               to be symmetric (evens out the colormapeditor step)
% The result (sampled at 64 color points) is:

hue_ranges = [ -5 5 6 10 11 14 15 27 28 35 36 50 52 59 ];

hue_vals = [];
for j = 1:(length(hue_ranges)/2)
    hue_vals = [hue_vals  linspace( hue_ranges( (2*j) -1 ), hue_ranges(2*j), 10 ) ];
end
hue_vals = mod(hue_vals/64,1)';

hue_vals = circshift(hue_vals,15);% i want colors to still start on blue

map = [ hue_vals ones(size(hue_vals)) ones(size(hue_vals))];
map = hsv2rgb(map);

map = interp1(linspace(0,1,size(map,1)), map, linspace(0,1,n), 'linear');