function IV = CommonIV(cfg,ivA,ivB,varargin)
%COMMONIV Return intervals that are common across iv structs. Default
%settings return the intersection of the intervals in time. Output is the
%same regardless of the order of the inputs. This function is slow.
%
%   IV = COMMONIV(cfg,ivA,ivB,varargin)
%
%   ivA _____________________     ______    _________________    ____
%   ivB   _________   __________                __________             ____
%
%   IV    _________   _______                   __________
%
%       IV contains the time intervals where ivA and ivB intersect if
%       cfg.keepGaps = 1 and cfg.threshold = 1.
%
%       For more information about what this function is doing, open it and
%       read the main body of the code (contains ASCII drawings of steps)
%
%   INPUTS
%      cfg: config struct with fields controlling function behaviour (see
%           CONFIG OPTIONS below)
%      ivA, ivB, and varargins: iv data with intervals sorted in ascending order
%
%   OUTPUTS
%      IV:  iv struct containing intervals that are common across the
%           input iv structs, as defined by the user. The intervals are
%           non-overlapping.
% 
%   CONFIG OPTIONS
%      cfg.binsize = 1/2000; Internal timebase (same units as the input IVs)
%                    for deciding whether the intervals intersect. If the
%                    binsize is too large, merging can occur.
%      cfg.threshold = 1; What proportion of inputs need to agree on a 
%                    timespan for it to be kept as an interval? (this is a
%                    number between 0 and 1, such as 2/3 or 9/10)
%
%      cfg.keepGaps = 1; If 1, preserves regions where at least one of the
%                    inputs considered to be separate; if 0, doesn't. Note
%                    that merging can be performed with MergeIV if desired.
%
%      cfg.verbose = 1; If 1, print informative text to the command window;
%                    if 0, be silent.  
%
% **Implementation could be improved: currently uses a binning method that
%   was quick to code but is slow and suboptimal (any ideas? will this ever
%   be read by another human? complete rewrite?)
%
%   see also: UnionIV2, OverlapIV
%
% aacarey Feb 2016

cfg_def.binsize = 1/2000;
cfg_def.threshold = 1;
cfg_def.keepGaps = 1; % if 1, preserves areas that any of the intervals think should be separate
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

if cfg.verbose; tic; fprintf('%s: looking for time intersections in iv data...\n',mfun); end

nArgs = 2 + length(varargin);

% we could do a CheckIV here, but it's redundant since MergeIV below also
% does it, though the user won't know specifically which iv is bad...

%~~~ don't allow overlap within individual input intervals ~~~
cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.gap = 0;
ivA = MergeIV(cfg_temp,ivA);
ivB = MergeIV(cfg_temp,ivB);
for iVarg = 1:length(varargin)
    varargin{iVarg} = MergeIV(cfg_temp,varargin{iVarg});
end

%~~~ use UnionIV2 to collect the inputs into a single iv struct ~~~
cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.rmdoubles = 0;
% rmdoubles is zero because if there ARE copies across ivs generated
% differently, then we want to keep this because it's an agreement
IV_temp = UnionIV2(cfg_temp,ivA,ivB,varargin{:});

% suppose we have these intervals:
%                ____     __               ____________________   
%          _________________                      _______
%       __________________________       ___________________________
%                      <---------- tvec ---------->

%~~~ find the earliest and latest entries in all the input IVs ~~~
earliest = IV_temp.tstart(1); % starts are already sorted ascending
latest = max(IV_temp.tend); % ends are not necessarily sorted

%~~~ construct internal tvec ~~~
tvec = earliest:cfg.binsize:latest;

%~~~ bin interval agreement, this pairs with tvec so we can extract boundaries later ~~~
commonBins = zeros(1,length(tvec));

for iInterval = 1:length(IV_temp.tstart)
    withinInterval = find(tvec >= IV_temp.tstart(iInterval) & tvec <= IV_temp.tend(iInterval));
    commonBins(withinInterval) = commonBins(withinInterval) + 1;
end

% now we have somthing like this in commonBins, where the numbers 1-4
% indicate how many intervals contained the time bin in question:
%                ____     __                      _______
% 3        _____|    |___|  |              ______|       |_____
% 2     __|                 |_____       _|                    |____
% 1 ___|                          |_____|                           |_____
%                      <---------- tvec ---------->

%~~~ threshold commonBins ~~~
% convert threshold into a number compared to the total number of input ivs
threshold = nArgs*cfg.threshold;
commonBins(commonBins < threshold) = 0;

% so now maybe we have something like this:
%                ____     __                      _______
% 3        _____|    |___|  |              ______|       |_____
% 2       |                 |             |                    |
% 1 ______|                 |_____________|                    |__________
%                      <---------- tvec ---------->

if cfg.keepGaps % find valleys and pull them down
    % invert commonBins so that the valleys are now peaks that we can
    % detect with findpeaks
    
    commonBinsInverted = commonBins * -1;
    [~,minHere] = findpeaks(commonBinsInverted);
    commonBins(minHere) = 0;
% 4              ____     __                      _______
% 3        _____|    |   |  |              ______|       |_____
% 2       |          |   |  |             |                    |
% 1 ______|          |___|  |_____________|                    |__________
%                      <---------- tvec ---------->
end

%~~~ convert to tsd and keep everything above zero ~~~
commonTSD = tsd(tvec,commonBins);
cfg_temp = []; cfg_temp.verbose = 0; cfg.method = 'raw'; cfg_temp.threshold = 0; cfg_temp.operation = '>';
IV = TSDtoIV2(cfg_temp,commonTSD);

% if we wanted 2/3 to agree, then our intervals look like:
%          __________     __               ____________________
%                      <---------- tvec ---------->

%~~~ adjust boundaries ~~~
% the boundaries in IV should be derived from the boundaries in IV_temp
IV.tstart = nearval3(IV.tstart,IV_temp.tstart,0);
% binary search requires that these are in sorted order:
lookups = sort(IV_temp.tend);
IV.tend = nearval3(IV.tend,lookups,0);

if cfg.verbose; fprintf('%s: %d intervals in, %d intervals out\n',mfun,length(IV_temp.tstart),length(IV.tstart)); toc; end

% write history
IV.cfg.history = IV_temp.cfg.history; % keep the history of the predecessor
IV = History(IV,mfun,cfg);

end
