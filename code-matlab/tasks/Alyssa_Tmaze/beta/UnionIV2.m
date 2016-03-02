function IV = UnionIV2(cfg,ivA,ivB,varargin)
% UNIONIV2 Returns the set of all intervals without exact duplicates,
% though overlap is allowed
%
% function iv = UNIONIV2(cfg,ivA,ivB,varargin)
%
%  ivA     _______       ___      _________      ___________
%  ivB       ________         __                 ___________
%
%
%  IV      _______       ___  __  _________      ___________
%            ________ 
%          (^ note overlap)                      (^ one copy)
%
%
%    INPUTS
%      cfg: config struct with fields controlling function behaviour (see
%           CONFIG OPTIONS below)
%      ivA, ivB, and varargins: iv data with intervals sorted in ascending
%           order. Note: intervals may overlap and since some other
%           functions assume that this is not the case, you should consider
%           merging the intervals (MergeIV).
%
%    OUTPUTS
%      IV: iv struct containing all unique intervals present in the input
%          ivs, resorted based on start times (ascending). If cfg.verbose,
%          the function will report the total number of incoming intervals
%          and the number of unique outgoing intervals.
%
%    CONFIG OPTIONS
%      cfg.rmdoubles = 1; A true union does not contain multiple copies,
%          but you can disable this feature if you need it for something by
%          setting this ("remove doubles") to 0.
%      cfg.verbose = 1; If 1, print informative text to the command window;
%                      if 0, be silent.
%
%  see also: CommonIV, MergeIV
%
% MvdM 2014-08-28 initial version
% aacarey edit Feb 2015, takes variable number of IV inputs, rejects duplicates

% set config defaults
cfg_def.rmdoubles = 1;
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

nArg = 2 + length(varargin);

% collect all IVs into a single variable for ease of looping
allIVs = cell(1,nArg);
allIVs{1} = ivA; allIVs{2} = ivB;
if ~isempty(varargin)
    for iArg = 3:nArg % collect the remaining IVs
        allIVs(iArg) = varargin(iArg-2);
    end
end

% check that the IVs are well formed
if any(cellfun(@(x) ~CheckIV(x), allIVs))
    error('At least one iv input is poorly formed')
end

% preallocate space so matlab doesn't orange at me (not sure if this is
% actually faster than concatenating them)
nIntervals = cellfun(@(x) length(x.tstart),allIVs);
totalIntervals = sum(nIntervals);
tstart = NaN(totalIntervals,1); tend = NaN(totalIntervals,1);

% collect all start and end times
next = 1; % tracks which index we're at
for iArg = 1:nArg
   % add the start and end times for this interval to the tstart and tend
   % collectors
   tstart(next:next+nIntervals(iArg)-1) = allIVs{iArg}.tstart;
   tend(next:next+nIntervals(iArg)-1) = allIVs{iArg}.tend;
   next = next + nIntervals(iArg);
end

% sort intervals and make output IV
[tstart,sort_idx] = sort(tstart,'ascend');
tend = tend(sort_idx);

IV = iv(tstart,tend);

% duplicate intervals are not allowed (exact matches are probably very rare)
if cfg.rmdoubles
    cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.mindur = 0; cfg_temp.maxdur = [];
    cfg_temp.rmdoubles = cfg.rmdoubles; % *** remove doubles
    IV = RemoveIV(cfg_temp,IV);
end

% talk to me
if cfg.verbose; fprintf('%s: %d intervals in, %d intervals out\n',mfun,totalIntervals,length(IV.tstart)); end

% write history
IV = History(IV,mfun,cfg);

end
