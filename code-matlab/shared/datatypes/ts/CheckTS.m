function pass_flag = CheckTS(ts_in,varargin)
%% CHECKTS Check TS for datatype violations
%	pass_flag = CheckTS(ts_in,varargin) verifies that input is ts and is well formed.
%
%	INPUTS:
%       ts_in: ts to be checked
%       varargins:
%           callername: name of invoking function (for more informative warning messages)
%
%	OUTPUTS:
%       pass_flag: 1 if all checks pass, 0 if otherwise
%
%	Checks performed:
%       -  .t field must exist (FAIL)
%       -  all contents of .t must be column vectors (FAIL)
%       -  is there a .label field? (WARNING)
%       -  is .label the same length as .t? (FAIL)
%       -  is there a .type field? (WARNING)
%       -  does the .type field correctly identify as ts? (WARNING)
%
%   see also ts, CheckTSD, CheckIV, CheckTC
%
% aacarey Nov 2015
% youkitan edit Dec 2016, reformat help, add function name to output

pass_flag = 1;

in_mfun = '';
if ~isempty(varargin) && ischar(varargin{1})
    in_mfun = [' in ',varargin{1}];
end

if ~isstruct(ts_in)
    pass_flag = 0;
    fprintf('FAIL%s by CheckTS: input must be a ts datatype.\n',in_mfun);
else
    % check for .t field
    if ~isfield(ts_in,'t')
        pass_flag = 0;
        fprintf('FAIL%s by CheckTS: input ts must contain t field.\n',in_mfun);
        
    else
        % check that all t are column vectors
        for iT = 1:length(ts_in.t)
            if ~iscolumn(ts_in.t{iT})
                pass_flag = 0;
                fprintf('FAIL%s by CheckTS: all contents of ts_in.t must be column vectors.\n',in_mfun);
            end
        end
    end
    
    % check for label field
    if ~isfield(ts_in,'label')
        %pass_flag = 0;
        fprintf('WARNING%s by CheckTS: ts_in lacks a .label field.\n',in_mfun);
    else
        % check that label is the same length as t
        if length(ts_in.label) ~= length(ts_in.t)
            pass_flag = 0;
            fprintf('FAIL%s by CheckTS: ts_in.t and ts_in.label must have the same length.\n',in_mfun);
        end
    end

    % check for type field
    if pass_flag && ~isfield(ts_in,'type')
        %pass_flag = 0;
        fprintf('WARNING%s by CheckTS: input appears to be ts datatype but lacks the .type field.\n',in_mfun);
    elseif pass_flag && isfield(ts_in,'type') && ~strcmp(ts_in.type,'ts')
        %pass_flag = 0;
        fprintf('WARNING%s by CheckTS: input appears to be ts datatype but is identified as %s in the .type field.\n',in_mfun,ts_in.type);
    end
    
end

end