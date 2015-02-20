function [left,right] = GetMatchedTrials2(cfg_in,metadata,ExpKeys)
%GetMatchedTrials2 Return equal subsets of trials for L & R
%
% First, this function removes bad trials from the trial sets.
%
% Next, it checks the number of remaining trials.
% If there's an equal number of L and R trials (after bad trials are
% considered), it returns all of the good trials.
% If one set outnumbers the other, it returns a random subset of trials from
% this set, and all of the trials from the lower-number set.
% *** the subset it returns is randomly generated, but repeated runs will
% always generate the same subset.
% 
% ACarey Feb 2015

cfg_def = [];

[~] = ProcessConfig2(cfg_def,cfg_in);

%% pull out left and right trials from metadata

left = metadata.taskvars.trial_iv_L;
right = metadata.taskvars.trial_iv_R;

%% account for bad (discarded) trials:

if ~isempty(ExpKeys.badTrials)

    Lkey = find(strcmp('L',metadata.taskvars.sequence));
    Rkey = find(strcmp('R',metadata.taskvars.sequence));
    
    % the above "keys" contain the original indices of the left and right
    % trials with respect to the whole trial sequence. They can later be
    % used to match the discarded trial's index (in the whole trial set)
    % with the L- or R-specific index in the L- or R-only sets. For example, 
    % if trial 18 is a discarded left trial according to the whole set
    % of trials, the number 18 will be present in the Lkey. The index that
    % 18 sits at in Lkey is the index of the L intervals we need to remove.
    
    Lkill = []; % to hold indices of bad L trials
    Rkill = []; %to hold indices of bad R trials
    
    for iDiscard = 1:length(ExpKeys.badTrials)
        %disp(ExpKeys.badTrials(iDiscard))
        if strcmp('L',metadata.taskvars.sequence(ExpKeys.badTrials(iDiscard)))
            Lkill = find(Lkey == ExpKeys.badTrials(iDiscard));
            disp('Discarding L trial (bad).')
        elseif strcmp('R',metadata.taskvars.sequence(ExpKeys.badTrials(iDiscard)))
            Rkill = find(Rkey == ExpKeys.badTrials(iDiscard));
            disp('Discarding R trial (bad).')
        end
    end
    
    % get rid of the trials
    if ~isempty(Lkill)
        left.tstart(Lkill) = []; left.tend(Lkill) = [];
    end
    
    if ~isempty(Rkill)
        right.tstart(Rkill) = []; right.tend(Rkill) = [];
    end
end

%% Return trial iv sets

% get trial counts
nLtrials = length(left.tstart);
nRtrials = length(right.tstart);

% generate subsets
if nLtrials == nRtrials
    fprintf('Equal good trials: nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    return;
    
elseif nLtrials > nRtrials
    fprintf('Good trials: nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    left = predsubset(left,nRtrials);
    
else % nRtrials > nLtrials
    fprintf('Good trials: nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    right = predsubset(right,nLtrials);
end
      
end

