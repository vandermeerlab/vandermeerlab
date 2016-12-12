function pass_flag = CheckIV(iv_in,varargin)
%% CHECKIV Check IV for datatype violations
%	pass_flag = CheckIV(ts_in,varargin) verifies that input is iv and is well formed.
%
%	INPUTS:
%       iv_in: iv to be checked
%       varargins:
%           callername: name of invoking function (for more informative warning messages)
%
%	OUTPUTS:
%       pass_flag: 1 if all checks pass, 0 if otherwise
%
%	Checks performed:
%       -  .tstart and .tend fields must exist (FAIL)
%       -  are either .tstart or .tend empty? (WARNING)
%       -  both .tstart and .tend must be column vectors (FAIL)
%       -  are start-end pairs unidirectional? (WARNING)
%
%	see also iv, CheckTS, CheckTSD, CheckTC
%
% MvdM 2014-11-12
% aacarey edit Sept 2015, additional checks
% youkitan edit Sept 2016, additional checks
% youkitan edit Dec 2016, reformat help, add function name to output

pass_flag = 1;

in_mfun = '';
if ~isempty(varargin) && ischar(varargin{1})
    in_mfun = [' in ',varargin{1}];
end

if isstruct(iv_in)
    if ~isfield(iv_in,'tstart') || ~isfield(iv_in,'tend')
        pass_flag = 0;
        fprintf('FAIL%s by CheckIV: input iv must contain tstart and tend fields.\n',in_mfun);
    elseif isempty(iv_in.tstart) || isempty(iv_in.tend)
        fprintf('WARNING%s by CheckIV: input iv is empty.\n',in_mfun);
    elseif ~iscolumn(iv_in.tstart) || ~iscolumn(iv_in.tend)
        pass_flag = 0;
        fprintf('FAIL%s by CheckIV: tstart and tend must be column vectors.\n',in_mfun);
    elseif any(iv_in.tstart > iv_in.tend)
        fprintf('WARNING%s by CheckIV: start-end pairs not unidirectional.\n',in_mfun);
    end
else
    pass_flag = 0;
    fprintf('FAIL%s by CheckIV: input must be an iv data type.\n',in_mfun);
end
   

