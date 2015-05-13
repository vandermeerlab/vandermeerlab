function pass_flag = CheckIV(iv_in)
% function pass_flag = CheckIV(iv_in)
%
% checks if iv is well formed
%
% INPUTS:
% iv_in: iv to be checked
%
% OUTPUTS:
% pass_flag: 1 if all checks pass, 0 if otherwise
%
% MvdM 2014-11-12

pass_flag = 1;

if ~iscolumn(iv_in.tstart) || ~iscolumn(iv_in.tend)
   pass_flag = 0;
   fprintf('FAIL: t_start and/or t_end are not column vectors.\n');
end
