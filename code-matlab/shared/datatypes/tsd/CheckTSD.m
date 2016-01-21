function pass_flag = CheckTSD(tsd_in)
% function pass_flag = CheckTSD(tsd_in)
%
% checks if TSD is well formed
%
% INPUTS:
% tsd_in: tsd to be checked
%
% OUTPUTS:
% pass_flag: 1 if all checks pass, 0 if otherwise
%
% MvdM 2014-06-24

pass_flag = 1;

if ~isstruct(tsd_in)
   pass_flag = 0;
   fprintf('FAIL: data is not a struct.\n');
   return;
end

if iscell(tsd_in.data)
   pass_flag = 0;
   fprintf('FAIL: data is a cell.\n');
   return;
end

nSignals = size(tsd_in.data,1);
nSamples_data = size(tsd_in.data,2);

if nSignals > nSamples_data
   pass_flag = 0;
   fprintf('WARNING: more signals (%d) than samples (%d) in data.\n',nSignals,nSamples_data);
end

nSamples_tvec = length(tsd_in.tvec);
if nSamples_data ~= nSamples_tvec
    pass_flag = 0;
    fprintf('FAIL: samples in data (%d) does not match samples in tvec.\n',nSamples_data,nSamples_tvec);
end