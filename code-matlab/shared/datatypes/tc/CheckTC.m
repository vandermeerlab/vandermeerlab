function pass_flag = CheckTC(tc_in,varargin)
%CHECKTS Verify that input is ts and is well formed
%
%    INPUTS:
%       tc_in: tc to be checked
% 
%       varargins:
%       callername: name of invoking function (for more informative warning
%       messages)
%
%       OUTPUTS:
%       pass_flag: 1 if all checks pass, 0 if otherwise
%
%      Checks performed:
%        -  .tc field must exist (FAIL)
%        -  binsize must agree across variables (FAIL)
%        -  number of cells must agree across variables (FAIL)
%
% youkitan 2016-11-29 initial

pass_flag = 1;

in_mfun = '';
if ~isempty(varargin) && ischar(varargin{1})
    in_mfun = [' in ',varargin{1}];
end

if ~isstruct(tc_in)
    pass_flag = 0;
    fprintf('FAIL%s: input must be a tc struct.\n',in_mfun);
else
    % checks for 1D case
    if isfield(tc_in,'tc')
        if size(tc_in.tc,2) ~= size(tc_in.occ_hist,2)
            pass_flag = 0;
            fprintf('FAIL%s: number of bins differ between tc and occ_hist.\n',in_mfun);
        end
        
    % checks for 2D case    
    elseif isfield(tc_in,'tc2D')
        fprintf('WARNING: Checks for 2D tuning curves not implemented yet!!!')
    
    else %no tuning variable
        pass_flag = 0;
        fprintf('FAIL%s: input tc must contain a tc field.\n',in_mfun);
    end
   
end

end

