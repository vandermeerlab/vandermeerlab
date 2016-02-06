function [left,right,left_indices,right_indices] = GetMatchedTrials(cfg_in,metadata,ExpKeys)
% function [left,right] = GetMatchedTrials(cfg_in,metadata,ExpKeys)
%
% left_indices and right_indices are the trial numbers for left and right
% trials in the full original sequence (which includes bad trials). 
%
% MvdM
% AC edit, handles bad trials

cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

%% account for bad (discarded) trials:

if ~isempty(ExpKeys.badTrials)
   
    sequence = metadata.taskvars.sequence;
    sequence(ExpKeys.badTrials) = [];
 
    trial_numbers = 1:length(metadata.taskvars.trial_iv.tstart);
    trial_numbers(ExpKeys.badTrials) = [];
    
    L_trial_numbers = trial_numbers(strcmp('L',sequence));
    
    R_trial_numbers = trial_numbers(strcmp('R',sequence));
    
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
    if cfg.verbose; fprintf('%s: Equal good trials: nLtrials %d, nRTrials %d...\n',mfun,nLtrials,nRtrials); end
    left = metadata.taskvars.trial_iv_L;
    right = metadata.taskvars.trial_iv_R;
    left_indices = L_trial_numbers;
    right_indices = R_trial_numbers;
    return;

elseif nLtrials > nRtrials
    if cfg.verbose; fprintf('%s: Good trial counts: nLtrials %d, nRTrials %d...\n',mfun,nLtrials,nRtrials); end
    match = 'R';
     
    % these contain a list of the pool of trial numbers for each trial type
    min_trial_no = R_trial_numbers;
    
    max_trial_no = L_trial_numbers;
    
else
    if cfg.verbose; fprintf('%s: Good trial counts: nLtrials %d, nRTrials %d...\n',mfun,nLtrials,nRtrials); end
    match = 'L';
    
    % these contain a list of the pool of trial numbers for each trial type
    min_trial_no = L_trial_numbers;
    
    max_trial_no = R_trial_numbers;
end


% choose nearest trials

matched_trial = zeros(size(min_trial_no));

for iT = 1:length(min_trial_no)
    
    max_trial_idx = nearest_idx3(min_trial_no(iT),max_trial_no); % nearest_idx errors when matching last idx 
    matched_trial(iT) = max_trial_no(max_trial_idx);
    
    max_trial_no(max_trial_idx) = []; % remove matched trial from available pool to prevent duplicates

end

%
if length(unique(matched_trial)) < length(matched_trial) % things somehow got messed up
    warning('*** Could not determine matched trial sequence.');
end

matched_trial = sort(matched_trial); % ensures they are temporally ordered

% return subsets
switch match
    case 'L'
        left.tstart = metadata.taskvars.trial_iv.tstart(L_trial_numbers);
        left.tend = metadata.taskvars.trial_iv.tend(L_trial_numbers);
        left = iv(left.tstart,left.tend);
       
        right.tstart = metadata.taskvars.trial_iv.tstart(matched_trial);
        right.tend = metadata.taskvars.trial_iv.tend(matched_trial);
        right = iv(right.tstart,right.tend);
        
        left_indices = L_trial_numbers;
        right_indices = matched_trial;
        
    case 'R'
        left.tstart = metadata.taskvars.trial_iv.tstart(matched_trial);
        left.tend = metadata.taskvars.trial_iv.tend(matched_trial);
        left = iv(left.tstart,left.tend);
        
        right.tstart = metadata.taskvars.trial_iv.tstart(R_trial_numbers);
        right.tend = metadata.taskvars.trial_iv.tend(R_trial_numbers);
        right = iv(right.tstart,right.tend);
        
        left_indices = matched_trial;
        right_indices = R_trial_numbers;
end

end

