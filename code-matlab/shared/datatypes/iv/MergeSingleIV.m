function iv_out = MergeSingleIV(iv_in)
%MERGEIV Merge touching or overlapping intervals within a single iv struct
%   If an interval's end time is less than or equal to the next interval's 
%   start time, the two intervals are merged into a single interval.
%
%   iv_out = MergeSingleIV(iv_in)
%
% ACarey, March 2015

is_iv = CheckIV(iv_in);

if is_iv
    kill = [];
    for iInterval = 1:length(iv_in.tstart)-1;
        if iv_in.tstart(iInterval+1) <= iv_in.tend(iInterval)
            kill = [kill iInterval];
            iv_in.tstart(iInterval+1) = iv_in.tstart(iInterval);
        end 
    end
    iv_out = iv_in;
    if ~isempty(kill)
        iv_out.tstart(kill) = [];
        iv_out.tend(kill) = [];
    end
end


