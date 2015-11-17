function [iv_out,idx] = SelectIV2(cfg_in,iv_in,selectspec)
%SELECTIV2 Specify intervals to keep.
%   [iv_out,idx] = SELECTIV2(cfg_in,iv_in,selectspec)
%
%   INPUTS:
%         cfg: config struct with fields controlling function behavior
%       iv_in: iv struct
%  selectspec: selection specifics, either:
%         - [nx1] double: logical array or indices specifying which
%                 intervals to keep.
%         - string: string specifying which usr field to work with. If
%                 selectspec is a string, the config options cfg.operation,
%                 cfg.threshold, and cfg.str apply.

%
%   OUTPUTS
%      iv_out: iv struct with specified intervals selected and all
%              corresponding same-length usr trimmed accordingly
%
%   CFG OPTIONS
%       cfg.operation = '>='; How to perform numerical selection, see
%                     cfg.threshold.
%       cfg.threshold = 0; Set a numerical threshold for keeping intervals.
%                     This works on numerical usr contents, but can also be
%                     applied to strings as long as the first character is
%                     number-convertible:
%                     If your field contains strings and cfg.str is
%                     empty, SelectIV assumes that the first character is a
%                     number (i.e. a rating) and thresholds based on this 
%                     number. An example would be '1, very good', for which
%                     SelectIV considers the 1 only. 
%       cfg.str = ''; If your target usr field contents contain strings that 
%                     are NOT number-convertible, input the string you
%                     want to select by. If this is not empty, it overrides 
%                     numerical selection. Examples of non-number-convertible 
%                     strings might be 'good' or 'maybe' or 'poor'.
%       cfg.verbose = 1; Tell me how many intervals came in, and how many
%                     went out.
% aacarey Nov 2015
%
% see also restrict SelectIV RemoveIV restrict

cfg_def.operation = '>=';
cfg_def.threshold = 0;
cfg_def.str = ''; % if this is not empty, it overrides numerical selection
cfg_def.verbose = 1;

if ~CheckIV(iv_in)
    error('iv_in must be an iv data type.')
end

mfun = mfilename;

% parse cfg parameters
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% choose which thing to do
if islogical(selectspec) || isnumeric(selectspec)
    goto = 'there'; % -_-
elseif ischar(selectspec)
    goto = 'here'; % D:
else
    error('selectspec must be a logical array, numeric array of indices, or a string specifying a usr field name.')
end

switch goto
    case 'here'  
        % make sure usr exists
        if ~isfield(iv_in,'usr')
            error([mfun,': iv_in requires usr for this type of selection.'])
        end
        % check that the field actually exists and that it's the right length
        if ischar(selectspec) && ~isfield(iv_in.usr,selectspec)
            error([selectspec,' does not exist.'])
        elseif ischar(selectspec) && length(iv_in.usr.(selectspec)) ~= length(iv_in.tstart)
            error(['iv_in.usr.',selectspec,' must have the same dimensions as iv_in.tstart.'])
        end
        
        % if the field contains strings, get ratings in numerical form
        if isempty(cfg.str) && ~isnumeric(iv_in.usr.(selectspec)(1))
            type = 1;
            temp = nan(size(iv_in.usr.(selectspec)));
            for ii = 1:length(temp)
                temp(ii,1) = str2double(iv_in.usr.(selectspec){ii,1}(1)); % we assume the rating is the first character in the string
            end
        elseif isempty(cfg.str) && isnumeric(iv_in.usr.(selectspec))
            type = 1;
            temp = iv_in.usr.(selectspec);
        elseif ~isempty(cfg.str)
            type = 2;
            temp = iv_in.usr.(selectspec);
        end
        
        assignin('base','temp',temp)
        
        % do the thing
        switch type
            case 1
                switch cfg.operation
                    case '>'
                        keep =  temp > cfg.threshold;
                    case '>='
                        keep = temp >= cfg.threshold;
                    case '<'
                        keep = temp < cfg.threshold;
                    case '<='
                        keep = temp <= cfg.threshold;
                    case '='
                        keep = temp == cfg.threshold;
                    otherwise
                        error('Unrecognized cfg.operation')
                end
            case 2
                keep = nan(size(temp));
                for iStr = 1:length(temp)
                    keep(iStr) = strcmp(cfg.str,temp(iStr));
                end
        end
        
        keep = logical(keep);
        
    case 'there'
        keep = selectspec;
end

iv_out = iv_in;
iv_out.tstart = iv_out.tstart(keep);
iv_out.tend = iv_out.tend(keep);

% also select data from other same-length usr fields
if isfield(iv_out,'usr') && ~isempty(iv_out.usr)
    ivfields = fieldnames(iv_out.usr);
    for iField = 1:length(ivfields)
        iv_out.usr.(ivfields{iField}) = iv_out.usr.(ivfields{iField})(keep);
    end
end

% make idx output
if islogical(keep)
    idx = find(keep);
elseif isnumeric(keep)
    idx = keep; 
end

% talk to me
if cfg.verbose
    disp([mfun,': ',num2str(length(iv_in.tstart)),' intervals in, ',num2str(length(iv_out.tstart)),' intervals out.'])
end

% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_in.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_in.cfg.history.cfg,{cfg});

end % of function

