function [x,x_idx] = trim_histc(x,varargin)
% function [x,x_idx] = trim_histc(x,x_idx)
%
% trims output of histc() to merge last bin with second-to-last bin
%
% x_idx is optional
%
% WARNING: do not use with 2-D inputs!
%
% MvdM 2014-08-21

if isrow(x)
    rowflag = 1;
    x = x';
else
    rowflag = 0;
end

if nargin == 2
    x_idx = varargin{1};
    x_idx(x_idx == length(x)) = length(x)-1;
else
   x_idx = []; 
end

x(end-1) = x(end-1)+x(end);
x = x(1:end-1);

if rowflag
   x = x'; 
end


