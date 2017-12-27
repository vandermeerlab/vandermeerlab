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
% cfdu@uwaterloo.ca Nov 2017 implemented to remove binning process, speed 500x

cfg_def.threshold = 1;
cfg_def.keepGaps = 1; % if 1, preserves areas that any of the intervals think should be separate
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

if cfg.verbose; tic; fprintf('%s: looking for time intersections in iv data...\n',mfun); end

args = {ivA, ivB, varargin{:}};
nArgs = length(args);

% we could do an isIV here, but it's redundant since MergeIV below also
% does it, though the user won't know specifically which iv is bad...

%~~~ don't allow overlap within individual input intervals ~~~
cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.gap = 0;
for iArg = 1:nArgs
    args{iArg} = MergeIV(cfg_temp,args{iArg});
end

% suppose we have these intervals:
%                ____     __               ____________________   
%          _________________                      _______
%       __________________________       ___________________________
%                      <---------- tvec ---------->

%                ____     __                      _______
% 3        _____|    |___|  |              ______|       |_____
% 2     __|                 |_____       _|                    |____
% 1 ___|                          |_____|                           |_____
%                      <---------- tvec ---------->

num_intervals = 0;
for iArg = 1:nArgs
    num_intervals = num_intervals + length(args{iArg}.tstart);
end
hist = nan(num_intervals, 2);
curr_index = ones(nArgs, 1);
curr_points = nan(nArgs, 1);
for iArg = 1:nArgs
    curr_points(iArg) = args{iArg}.tstart(1);
end

for hist_idx = 1:2*num_intervals
    [pos, min_idx] = min(curr_points);
    is_start = mod(curr_index(min_idx), 2);
    hist(hist_idx, 1) = pos;
    hist(hist_idx, 2) = (-1) ^ ~is_start;
    curr_index(min_idx) = curr_index(min_idx) + 1;
    next_arg_index = ceil(curr_index(min_idx) / 2);
    if is_start
        curr_points(min_idx) = args{min_idx}.tend(next_arg_index);
    else
        if next_arg_index <= length(args{min_idx}.tstart)
            curr_points(min_idx) = args{min_idx}.tstart(next_arg_index);
        else
            curr_points(min_idx) = inf;
        end
    end
end

hist(:, 2) = cumsum(hist(:,2));
hist = hist(logical([diff(hist(:, 1)); size(hist, 1)]), :);


%if cfg.keepGaps % find valleys and pull them down
% 4              ____     __                      _______
% 3        _____|    |   |  |              ______|       |_____
% 2       |          |   |  |             |                    |
% 1 ______|          |___|  |_____________|                    |__________
%                      <---------- tvec ---------->
if cfg.keepGaps
    temp = diff(hist(:, 2));
    hist([0; temp < 0] & [temp > 0; 0], 2) = 0;
end

% if we wanted 2/3 to agree, then our intervals look like:
%          __________     __               ____________________
%                      <---------- tvec ---------->
threshold = nArgs*cfg.threshold;
start_idx = find(hist(:, 2) >= threshold);
end_idx = start_idx + 1;
IV = iv(hist(start_idx), hist(end_idx));

if cfg.verbose; fprintf('%s: %d intervals in, %d intervals out\n',mfun,num_intervals,length(IV.tstart)); toc; end

% write history
IV = History(IV,mfun,cfg);

end
