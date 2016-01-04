function [ts_out,idx] = SelectTS(cfg_in,ts_in,selectspec)
%SELECTTS Specify ts data to keep. (Also performs reordering).
% This function may be of particular interest for ordering and selecting
% spiketrains (S output from LoadSpikes) but is applicable to any ts data.
%
%   ts_out = SELECTTS(cfg,ts_in,selectspec)
%
%   INPUTS:
%         cfg: config struct with fields controlling function behavior
%       ts_in: ts struct
%  selectspec: selection specifics, either:
%         - [nx1] double: logical array or indices specifying which
%                 intervals to keep.
%         - string: string specifying which usr field to work with. If
%                 selectspec is a string, the config options cfg.operation,
%                 cfg.threshold, and cfg.str apply.
%
%   OUTPUTS
%      ts_out - ts data ordered or selected according to selectspec
%
%   CONFIG OPTIONS
%       cfg.operation = '>='; How to perform numerical selection, see
%                     cfg.threshold.
%       cfg.threshold = 0; Set a numerical threshold for keeping ts.
%                     This works on numerical usr contents, but can also be
%                     applied to strings as long as the first character is
%                     number-convertible:
%                     If your field contains strings and cfg.str is
%                     empty, SelectTS assumes that the first character is a
%                     number (i.e. a rating) and thresholds based on this 
%                     number. An example would be '1, very good', for which
%                     SelectTS considers the 1 only. 
%       cfg.str = ''; If your target usr field contents contain strings that 
%                     are NOT number-convertible, input the string you
%                     want to select by. If this is not empty, it overrides 
%                     numerical selection. Examples of non-number-convertible 
%                     strings might be 'good' or 'maybe' or 'poor'.
%      cfg.verbose = 1; % if 1, tell me how many ts went in and how many 
%                went out; if 0, don't
%
% aacarey Oct 2015, edit Nov 2015 

%%


mfun = mfilename; 
if ~CheckTS(ts_in)
    error('ts_in is either not a ts datatype or is poorly formed.')
end

% Parse cfg parameters
cfg_def.operation = '>=';
cfg_def.threshold = 0;
cfg_def.str = ''; % if this is not empty, it overrides numerical selection
cfg_def.verbose = 1;

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
        if ~isfield(ts_in,'usr')
            error([mfun,': ts_in requires usr for this type of selection.'])
        end
        % check that the field actually exists and that it's the right length
        if ischar(selectspec) && ~isfield(ts_in.usr,selectspec)
            error([selectspec,' does not exist.'])
        elseif ischar(selectspec) && length(ts_in.usr.(selectspec)) ~= length(ts_in.t)
            error(['ts_in.usr.',selectspec,' must have the same dimensions as ts_in.t.'])
        end
        
        % if the field contains strings, get ratings in numerical form
        if isempty(cfg.str) && ~isnumeric(ts_in.usr.(selectspec)(1))
            type = 1;
            temp = nan(size(ts_in.usr.(selectspec)));
            for ii = 1:length(temp)
                temp(ii,1) = str2double(ts_in.usr.(selectspec){ii,1}(1)); % we assume the rating is the first character in the string
            end
        elseif isempty(cfg.str) && isnumeric(ts_in.usr.(selectspec))
            type = 1;
            temp = ts_in.usr.(selectspec);
        elseif ~isempty(cfg.str)
            type = 2;
            temp = ts_in.usr.(selectspec);
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

% select
ts_out = ts_in;
ts_out.t = ts_out.t(keep);
ts_out.label = ts_out.label(keep);

% also select data from other same-length usr fields
if isfield(ts_out,'usr') && ~isempty(ts_out.usr)
    TSfields = fieldnames(ts_out.usr);
    for iField = 1:length(TSfields)
        ts_out.usr.(TSfields{iField}) = ts_out.usr.(TSfields{iField})(keep);
    end
end

% make idx output
if islogical(keep)
    idx = find(keep);
elseif isnumeric(keep)
    idx = keep; 
end

if cfg.verbose % talk to me
    disp([mfun,': ',num2str(length(ts_in.t)),' ts in, ',num2str(length(ts_out.t)),' ts out'])
end

% keep a record of cfg history

if isfield(ts_in,'cfg')
    ts_out.cfg = ts_in.cfg;
    ts_out.cfg.history.mfun = cat(1,ts_in.cfg.history.mfun,mfilename);
    ts_out.cfg.history.cfg = cat(1,ts_in.cfg.history.cfg,{cfg});
else
    ts_out.cfg.history.mfun = mfun;
    ts_out.cfg.history.cfg = {cfg};
end

end

