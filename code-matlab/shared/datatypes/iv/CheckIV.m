function pass_flag = CheckIV(iv_in,varargin)
% function pass_flag = CheckIV(iv_in)
%
% checks if iv is well formed
%
% INPUTS:
% iv_in: iv to be checked
% 
% varargins:
%  callername: name of invoking function (for more informative warning
%  messages)
%
% OUTPUTS:
% pass_flag: 1 if all checks pass, 0 if otherwise
%
% MvdM 2014-11-12
% aacarey edit Sept 2015, additional checks

pass_flag = 1;

in_mfun = '';
if ~isempty(varargin) && ischar(varargin{1})
    in_mfun = [' in ',varargin{1}];
end

if isstruct(iv_in)
    if ~isfield(iv_in,'tstart') || ~isfield(iv_in,'tend')
        pass_flag = 0;
        fprintf('FAIL%s: input iv must contain tstart and tend fields.\n',in_mfun);
    elseif isempty(iv_in.tstart) || isempty(iv_in.tend)
        fprintf('WARNING%s: input iv is empty.\n',in_mfun);
    elseif ~iscolumn(iv_in.tstart) || ~iscolumn(iv_in.tend)
        pass_flag = 0;
        fprintf('FAIL%s: tstart and tend must be column vectors.\n',in_mfun);
    end
else
    pass_flag = 0;
    fprintf('FAIL%s: input must be an iv data type.\n',in_mfun);
end
   

