function pass_flag = CheckTSD(tsd_in, varargin)
%% CHECKTSD Check TSD for datatype violations
%   pass_flag = CheckTSD(tsd_in,varargin) verifies that input is TSD and well formed.
%
%	INPUTS:
%       tsd_in: tsd to be checked
%       varargins:
%           callername: name of invoking function (for more informative warning messages)
%
%	OUTPUTS:
%       pass_flag: 1 if all checks pass, 0 if otherwise
%
%   Checks performed:
%       -  .units field exists (WARNING)
%       -  .data field cannot be a cell (FAIL)
%       -  is the number of signals greater than the number of samples? (WARNING)
%       -  number of samples in .data and .tvec must be equal (FAIL)
%
%   see also tsd, CheckTS, CheckIV, CheckTC
%
% MvdM 2014-06-24
% youkitan edit Dec 2016, reformat help, add function name to output
% youkitan edit Feb 2017, included check for units

pass_flag = true;

in_mfun = '';
if ~isempty(varargin) && ischar(varargin{1})
    in_mfun = [' in ',varargin{1}];
end

if ~isstruct(tsd_in)
   pass_flag = false;
   fprintf('FAIL%s by CheckTSD: data is not a struct.\n', in_mfun);
   return;
end

if ~isfield(tsd_in,'units')
    fprintf('WARNING%s by CheckTSD: Missing units. TSD was generated with an outdated version of the constructor.\n', in_mfun);
end

if iscell(tsd_in.data)
   pass_flag = 0;
   fprintf('FAIL%s by CheckTSD: data is a cell.\n', in_mfun);
   return;
end

nSignals = size(tsd_in.data, 1);
nSamples_data = size(tsd_in.data, 2);

if nSignals > nSamples_data
   fprintf('WARNING%s by CheckTSD: more signals (%d) than samples (%d) in data.\n', in_mfun, nSignals, nSamples_data);
end

nSamples_tvec = length(tsd_in.tvec);
if nSamples_data ~= nSamples_tvec
    pass_flag = false;
    fprintf('FAIL in %s by CheckTSD: samples in data (%d) does not match samples in tvec.\n', in_mfun, nSamples_data, nSamples_tvec);
end