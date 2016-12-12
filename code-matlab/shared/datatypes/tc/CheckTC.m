function pass_flag = CheckTC(tc_in,varargin)
%% CHECKTC Check TC for datatype violations
%	pass_flag = CheckTC(tc_in,varargin) verifies that input is tc and is well formed.
%
%	INPUTS:
%       tc_in: tc to be checked
%       varargins:
%           callername: name of invoking function (for more informative warning messages)
%
%	OUTPUTS:
%       pass_flag: 1 if all checks pass, 0 if otherwise
%
%	Checks performed:
%       -  .tc or .tc2D field must exist (FAIL)
%       -  .tc must be an nCells x nBins matrix (FAIL)
%       -  .tc2D must be an nCells x nXBins x nYBins matrix (FAIL)
%       -  nBins must agree across variables (FAIL)
%       -  number of cells must agree across variables (FAIL)
%
%	see also: tc, TuningCurves, CheckTSD, CheckTS, CheckIV
%
% youkitan 2016-11-29 initial
% youkitan edit Dec 2016, added some tc2D functionality, reformat help

pass_flag = 1;

in_mfun = '';
if ~isempty(varargin) && ischar(varargin{1})
    in_mfun = [' in ',varargin{1}];
end

if ~isstruct(tc_in)
    pass_flag = 0;
    fprintf('FAIL%s by CheckTC: input must be a tc struct.\n',in_mfun);
else
    % basic checks
    if isfield(tc_in,'type')
        if ~isequal(tc_in.type,'tc')
            fprintf('FAIL%s by CheckTC: input must be a tc type.\n',in_mfun);
            pass_flag = 0;
            
        % checks for 1D case
        elseif isfield(tc_in,'tc')
            if length(size(tc_in.tc)) ~= 2
                fprintf('FAIL%s by CheckTC: tc variable should be an nCells x nBins matrix.\n',in_mfun);
                pass_flag = 0;

            elseif size(tc_in.tc,2) ~= size(tc_in.occ_hist,2)
                pass_flag = 0;
                fprintf('FAIL%s by CheckTC: number of bins differ between tc and occ_hist.\n',in_mfun);
            end

        % checks for 2D case    
        elseif isfield(tc_in,'tc2D')
            if length(size(tc_in.tc2D)) ~= 3
                fprintf('FAIL%s by CheckTC: tc variable should be an nCells x nXBins x nYBins matrix.\n',in_mfun);
                pass_flag = 0;
            end
            fprintf('WARNING by CheckTC: Checks for 2D tuning curves not implemented yet!!!')

        else %no tuning variable
            pass_flag = 0;
            fprintf('FAIL%s by CheckTC: input tc must contain a tc field.\n',in_mfun);
        end
    else
       fprintf('FAIL%s by CheckTC: input has no type parameter.\n',in_mfun)
       pass_flag = 0;
    end
   
end

end