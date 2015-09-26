function struct_in = structfunS(function_handle,struct_in,varargin)
% struct_out = structfunS(function_handle,struct_in,fieldnames)
%
% applies function to specific fields in struct, modifying them in place
%
% INPUTS:
% function_handle: [function handle] function to apply
% struct_in: [struct] input struct to run function on (can be 1-D struct array)
% (optional) fieldnames: [1 x nFields cell array] fields to apply function
% to; if not specified, apply to all
%
% OUTPUTS:
% struct_out: [struct] modified struct_in
%
% EXAMPLE USAGE:
% s.a = 1; s.b = 2; s.c = 3;
% fh = @(x) 1 + sqrt(x);
% s = structfunS(fh,s,{'b','c'});
%
% See also: structfun
%
% MvdM 2015-09-26 initial version

if ~isstruct(struct_in)
    error('Input is not a struct.');
end

if sum(size(struct_in) > 1) > 1 % count number of non-singleton dimensions
   error('Not yet implemented for >1D arrays.'); 
end

fn = fieldnames(struct_in(1));
if nargin == 3 % use specified fields only
   fn_todo = varargin{1};
   
   if length(fn_todo) < length(intersect(fn,fn_todo))
       error('One or more specified field names not found.');
   end 
else % use all fields
   fn_todo = fn; 
end

% loop over struct array elements
nel = length(struct_in);
for ii = 1:nel
    for iF = 1:length(fn_todo)
    
        struct_in(ii).(fn_todo{iF}) = function_handle(struct_in(ii).(fn_todo{iF}));
        
    end
end

