function iv_out = MergeSingleIV(cfg_in,iv_in)
%MERGESINGLEIV Merge touching or overlapping intervals within a single iv struct
%
%   iv_out = MergeSingleIV(cfg_in,iv_in)
%
%   If an interval's end time is less than or equal to the next interval's 
%   start time, the two intervals are merged into a single interval.
%
%   tstart(1) ______________________tend(1)
%                 tstart(2)______________________tend(2)
%
%   INPUTS:
%       iv_in: iv struct with intervals in ascending order
%
%   CFG OPTIONS: (none yet)
%
% ACarey, March 2015

cfg_def = [];
% a possible future config option would be merging intervals if they are
% within a certain user-specified distance from each other
cfg = ProcessConfig2(cfg_def,cfg_in);

% make sure the input was created using the iv constructor
is_iv = CheckIV(iv_in);

if is_iv % proceed
    mfun = mfilename;
    
    % make sure the intervals are ordered 
    if ~issorted(iv_in.tstart);
        error('Intervals must be in ascending order')
    end
    
    % do the merging
    kill = zeros(size(iv_in.tstart)); 
    for iInterval = 1:length(iv_in.tstart)-1;
        if iv_in.tstart(iInterval+1) <= iv_in.tend(iInterval)
            kill(iInterval) = iInterval;
            iv_in.tstart(iInterval+1) = iv_in.tstart(iInterval);
        end 
    end
    
    iv_out = iv_in;
    kill = kill > 0; % kill will be 1 for intervals to remove, and zero for the ones to be kept
    iv_out.tstart(kill) = [];
    iv_out.tend(kill) = [];

    % housekeeping (keep a record of cfg history)
    iv_out.cfg.history.mfun = cat(1,iv_out.cfg.history.mfun,mfun);
    iv_out.cfg.history.cfg = cat(1,iv_out.cfg.history.cfg,{cfg});
end

% Elyot's version of the merging step -_-
%iv_out.tstart = [iv_in.tstart(1) ; iv_in.tstart(iv_in.tstart > circshift(iv_in.tend,[1,0]))];
%iv_out.tend = [iv_in.tend(iv_in.tend < circshift(iv_in.tstart,[-1,0]));iv_in.tend(end)];



