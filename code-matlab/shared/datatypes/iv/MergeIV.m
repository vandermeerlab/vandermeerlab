function iv_out = MergeIV2(cfg_in,iv_in)
%MERGEIV Merge touching, overlapping, or nearby intervals within an iv struct
%   iv_out = MERGEIV(cfg,iv_in)
% 
%   If an interval's end time is less than or equal to the next interval's 
%   start time, the two intervals are merged into a single interval.
%
%               tstart(1) ______________________tend(1)
%  iv_in
%                             tstart(2)______________________tend(2)        
%
%             
%  iv_out       tstart(1) ___________________________________tend(1)        
%
%  If intervals are less than a user-defined distance away from each other, 
%  the intervals are merged (see cfg.gap).
%
%               tstart(1) ___________tend(1)
%
%   iv_in                         |<---cfg.gap--->|
%
%                                      tstart(2)_____________tend(2)        
%
%             
%   iv_out      tstart(1) ___________________________________tend(1)        
%
%   INPUTS:
%         cfg: config struct with fields controlling function behavior
%       iv_in: iv struct with intervals in ascending order
%
%   CFG OPTIONS: 
%       cfg.gap = 0; Merge events closer than this
%       cfg.verbose = 1; Tell me how many intervals came in, and how many
%                        went out.
% 
% aacarey oct 2015
%
%   see also UnionIV, IntersectIV, ResizeIV

%%
% set cfg defaults
cfg_def.gap = 0;
cfg_def.verbose = 1;

if ~CheckIV(iv_in);
    error('iv_in must be an iv data type.')
end

% make sure the intervals are ordered
if ~issorted(iv_in.tstart);
    error('Intervals must be in ascending order.')
end

mfun = mfilename;

% parse cfg parameters
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% add 1/2 cfg.gap to the intervals
cfg_temp.d = [-cfg.gap/2 cfg.gap/2];
cfg_temp.allowOverlap = 1;
iv_temp = ResizeIV(cfg_temp,iv_in);

% merge the resulting iv
remove = zeros(size(iv_temp.tstart));
for iInterval = 1:length(iv_temp.tstart)-1;
    if iv_temp.tstart(iInterval+1) <= iv_temp.tend(iInterval)
        remove(iInterval) = iInterval;
        iv_temp.tstart(iInterval+1) = iv_temp.tstart(iInterval);
    end
end

remove = remove > 0; % kill will be 1 for intervals to remove, and zero for the ones to be kept
iv_temp.tstart(remove) = [];
iv_temp.tend(remove) = [];

% remove excess from outer flanks (this happens because we're using ResizeIV
% to grow the inner flanks)
cfg_temp =[];
cfg_temp.d = [cfg.gap/2 -cfg.gap/2];
iv_out = ResizeIV(cfg_temp,iv_temp);

% talk to me
if cfg.verbose
    disp([mfun,': ',num2str(length(iv_in.tstart)),' intervals in, ',num2str(length(iv_out.tstart)),' intervals out.'])
end

% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_in.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_in.cfg.history.cfg,{cfg});
 
end

