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
% aacarey edit Sept 2015, additional checks

pass_flag = 1;

if isstruct(iv_in)
    if ~isfield(iv_in,'tstart') || ~isfield(iv_in,'tend')
        pass_flag = 0;
        fprintf('FAIL: input iv must contain tstart and tend fields.\n');
    end
    
    if ~iscolumn(iv_in.tstart) || ~iscolumn(iv_in.tend)
        pass_flag = 0;
        fprintf('FAIL: tstart and tend must be column vectors.\n');
    end
else
    pass_flag = 0;
    fprintf('FAIL: input must be an iv data type with tstart and tend fields.\n');
end
   

