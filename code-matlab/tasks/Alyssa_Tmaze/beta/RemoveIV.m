function iv_out = RemoveIV(cfg_in,iv_in)
%REMOVEIV Remove intervals that are too short or too long. Also checks for
%doubles with the option of removal. Corresponding usr data is also
%removed.
%   iv_out = REMOVEIV(cfg,iv_in)
%
%   INPUTS:
%         cfg: config struct with fields controlling function behavior
%       iv_in: iv struct 
%
%   OUTPUTS
%      iv_out: iv struct with specified intervals removed
%
%   CFG OPTIONS
%       cfg.mindur = 0; Minimum duration of intervals
%       cfg.maxdur = []; Maximum duration of intervals
%       cfg.rmdoubles = 0; 1 - remove doubles, 0 - don't
%       cfg.verbose = 1; If 1, tell me how many intervals came in, and how
%                        many went out. If 0, don't.   
%   Note: if you somehow got stuck in a paradox and have intervals that end 
%   before they begin, they will be removed.
%
% aacarey oct 2015
%
%   see also SelectIV
%%
% set cfg defaults
cfg_def.mindur = 0;
cfg_def.maxdur = [];
cfg_def.rmdoubles = 0;
cfg_def.verbose = 1;

if ~CheckIV(iv_in)
    error('iv_in must be an iv datatype.')
end

mfun = mfilename;

% parse cfg parameters
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

iv_temp = iv_in;
cfg_temp.verbose = 0; % prevent SelectIV from talking inside of RemoveIV

% check that there aren't any doubles
checkstart = find(diff(iv_temp.tstart)== 0);
checkend = find(diff(iv_temp.tend) == 0);
if any(checkend ~= checkstart) && cfg.verbose
    disp(['WARNING in ',mfun,': Some intervals have the same start or end times.'])
elseif ~isempty(checkstart) && cfg.rmdoubles
    if cfg.verbose
        disp([mfun,': ',num2str(length(checkstart>0)), ' doubles found, removing.'])
    end
    iv_temp.tstart(checkstart) = [];
    iv_temp.tend(checkstart) = [];
    
elseif ~isempty(checkstart) && ~cfg.rmdoubles && cfg.verbose
    disp(['WARNING in ',mfun,': ',num2str(length(checkstart>0)),' doubles found, but removal was not requested.'])
end

% removal based on duration of intervals
dur = iv_temp.tend - iv_temp.tstart;

if isempty(cfg.maxdur)
    keep = dur > cfg.mindur;
else
    keep = dur > cfg.mindur & dur < cfg.maxdur;
end

% notify if intervals end before they begin (or begin and end at same time)
if cfg.verbose
   oops = dur(dur <= 0); 
   if ~isempty(oops)
      disp(['OOPS in ',mfun,': the physics police have found and removed ',num2str(length(oops)),' intervals that end before or when they begin.']) 
   end
end

% make output
iv_out = SelectIV(cfg_temp,iv_temp,keep);
% iv_out = iv_temp;
% iv_out.tstart = iv_temp.tstart(keep);
% iv_out.tend = iv_temp.tend(keep);
% 
% % also remove data from other same-length usr fields
% if isfield(iv_in,'usr')
%     ivfields = fieldnames(iv_in.usr);
%     for iField = 1:length(ivfields)
%         iv_out.usr.(ivfields{iField}) = iv_out.usr.(ivfields{iField})(keep);
%     end
% end

% tell me how many
if cfg.verbose
    disp([mfun,': ',num2str(length(iv_in.tstart)),' intervals in, ',num2str(length(iv_out.tstart)),' intervals out.'])
end

% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_temp.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_temp.cfg.history.cfg,{cfg});

end

