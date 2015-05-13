function [left,right] = GetMatchedTrials(cfg_in,metadata)
% function [left,right] = GetMatchedTrials(cfg_in,metadata)

cfg_def = [];

cfg = ProcessConfig2(cfg_def,cfg_in);

nLtrials = length(metadata.taskvars.trial_iv_L.tstart);
nRtrials = length(metadata.taskvars.trial_iv_R.tstart);

match = '';
if nLtrials == nRtrials
    fprintf('Trial numbers equal (%d), trials returned as is.\n',nLtrials);
    left = metadata.taskvars.trial_iv_L;
    right = metadata.taskvars.trial_iv_R;
    return;
elseif nLtrials > nRtrials
    fprintf('nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    match = 'R';
    
    min_trial_iv = metadata.taskvars.trial_iv_R;
    min_trial_no = find(strcmp('R',metadata.taskvars.sequence));
    
    max_trial_iv = metadata.taskvars.trial_iv_L;
    max_trial_no = find(strcmp('L',metadata.taskvars.sequence));
else
    fprintf('nLtrials %d, nRTrials %d...\n',nLtrials,nRtrials);
    match = 'L';
    
    min_trial_iv = metadata.taskvars.trial_iv_L;
    min_trial_no = find(strcmp('L',metadata.taskvars.sequence));
    
    max_trial_iv = metadata.taskvars.trial_iv_R;
    max_trial_no = find(strcmp('R',metadata.taskvars.sequence));
end

%
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

%
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

