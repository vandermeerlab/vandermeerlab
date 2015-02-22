function [left,right] = GetMatchedTrials(cfg_in,metadata,ExpKeys)
% function [left,right] = GetMatchedTrials(cfg_in,metadata)


% MvdM
% AC edit, handles bad trials

cfg_def = [];

[~] = ProcessConfig2(cfg_def,cfg_in);


%% account for bad (discarded) trials:

if ~isempty(ExpKeys.badTrials)
    
    sequence = metadata.taskvars.sequence;
    sequence(ExpKeys.badTrials) = [];
    
    L_trial_numbers = find(strcmp('L',sequence));
    R_trial_numbers = find(strcmp('R',sequence));
    
else 
    
    sequence = metadata.taskvars.sequence;

    L_trial_numbers = find(strcmp('L',sequence));
    R_trial_numbers = find(strcmp('R',sequence));
    
end

%% Return trial iv sets

% get trial counts
nLtrials = length(L_trial_numbers);
nRtrials = length(R_trial_numbers);

% choose which set to downsize
if nLtrials == nRtrials
    fprintf('Equal good trials: nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    left = metadata.taskvars.trial_iv_L;
    right = metadata.taskvars.trial_iv_R;
    return;

elseif nLtrials > nRtrials
    fprintf('Good trial counts: nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    match = 'R';
     
    min_trial_no = R_trial_numbers;
    
    max_trial_no = L_trial_numbers; 
    
else
    fprintf('Good trial counts: nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    match = 'L';
    
    min_trial_no = L_trial_numbers;
    
    max_trial_no = R_trial_numbers;
end


% choose nearest trials

matched_trial = zeros(size(min_trial_no));

for iT = 1:length(min_trial_no)
    
    this_no = min_trial_no(iT);
    
    max_trial_idx = nearest_idx3(this_no,max_trial_no); % nearest_idx errors when matching last idx 
    matched_trial(iT) = max_trial_no(max_trial_idx);
    
    max_trial_no(max_trial_idx) = []; % remove matched trial from available pool to prevent duplicates
    
end

%
if length(unique(matched_trial)) < length(matched_trial) % things somehow got messed up
    warning('*** Could not determine matched trial sequence.');
end

% choose subsets
switch match
    case 'L'
        left = metadata.taskvars.trial_iv_L;
        
        right = metadata.taskvars.trial_iv;
        right.tstart = right.tstart(matched_trial); right.tend = right.tend(matched_trial);
        
    case 'R'
        left = metadata.taskvars.trial_iv;
        left.tstart = left.tstart(matched_trial); left.tend = left.tend(matched_trial);
        
        right = metadata.taskvars.trial_iv_R;
end

end

