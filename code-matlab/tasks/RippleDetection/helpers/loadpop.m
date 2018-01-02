function var_out = loadpop(filename,varargin)
%LOADPOP Load and rename a variable, and "pop it out"
%   Normally you can do var_out = load('fname','target_var_name');
%
% if the variable is called 'evt', then you have to access the information
% by doing var_out.evt.whatever
%
% this function pops evt out, so that you access it by just using
% var_out.whatever
%
% loadpop assumes that there is only a single variable stored inside of the
% file being loaded if target is unspecified
% if you specify the variable name to load in varargin, it will choose that one
%
% aacarey Nov 2015

if nargin == 1
    var_out = load(filename);
elseif nargin == 2
    var_out = load(filename,varargin{1});
else
    error('Unrecognized number of input arguments')
end

name = fieldnames(var_out); 

if length(name) > 1
  warning(['More than one variable detected in ',filename,' but target was unspecified'])
end

var_out = var_out.(name{1});

end

