function [metadata_out] = RemoveBadTrials(cfg_in,metadata_in,ExpKeys)
%REMOVEBADTRIALS send out identical copy of metadata, except that bad
%trials have been removed from the taskvars
%
%*** does not do anything with rest periods, so the outgoing sequence can't
%be used with the rest periods.
%
% Do not use GetMatchedTrials() on the output from this function. Use it on
% original metadata.
%
% ACarey Feb 2015. T-maze project

cfg_def = [];

[~] = ProcessConfig2(cfg_def,cfg_in);

if ~isempty(ExpKeys.badTrials)
   
    sequence = metadata_in.taskvars.sequence;
    sequence(ExpKeys.badTrials) = [];
 
    trial_numbers = 1:length(metadata_in.taskvars.trial_iv.tstart);
    trial_numbers(ExpKeys.badTrials) = [];
    
    L_trial_numbers = trial_numbers(strcmp('L',sequence));
    
    R_trial_numbers = trial_numbers(strcmp('R',sequence));
    
else 
    
    metadata_out = metadata_in;
    return
    
end

metadata_out = metadata_in;

metadata_out.taskvars.trial_iv.tstart =  metadata_in.taskvars.trial_iv.tstart(trial_numbers);
metadata_out.taskvars.trial_iv.tend =  metadata_in.taskvars.trial_iv.tend(trial_numbers);

metadata_out.taskvars.trial_iv_L.tstart =  metadata_in.taskvars.trial_iv.tstart(L_trial_numbers);
metadata_out.taskvars.trial_iv_L.tend =  metadata_in.taskvars.trial_iv.tend(L_trial_numbers);

metadata_out.taskvars.trial_iv_R.tstart =  metadata_in.taskvars.trial_iv.tstart(R_trial_numbers);
metadata_out.taskvars.trial_iv_R.tend =  metadata_in.taskvars.trial_iv.tend(R_trial_numbers);

metadata_out.taskvars.sequence = sequence;

end

