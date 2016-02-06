function metadata = TrimTrialTimes(cfg_in,metadata)
% function trial_iv = TrimTrialTimes(cfg,metadata)
%
%

cfg_def = [];
cfg_def.verbose = 0; % print output
cfg_def.mode = 'times'; % 'times', 'evt' (times.mat is specific to R042, evt mimics that)
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

switch cfg.mode
    
    case 'times'
        
        load(FindFile('*times.mat'));
        
        nT = length(metadata.taskvars.trial_iv.tstart);
        for iT = nT:-1:1
            
            start_dt(iT) = metadata.taskvars.trial_iv.tstart(iT)-run_start(iT);
            
            if start_dt(iT) < 0; % position-based trial starts before we have data
                metadata.taskvars.trial_iv.tstart(iT) = run_start(iT);
            end
            
            end_dt(iT) = metadata.taskvars.trial_iv.tend(iT)-run_end(iT);
            
            if start_dt(iT) > 0; % position-based trial ends after data cutoff in times
                metadata.taskvars.trial_iv.tend(iT) = run_end(iT);
            end
            
            left_idx = find(strcmp(metadata.taskvars.sequence,'L'));
            metadata.taskvars.trial_iv_L.tstart = metadata.taskvars.trial_iv.tstart(left_idx);
            metadata.taskvars.trial_iv_L.tend = metadata.taskvars.trial_iv.tend(left_idx);
            
            right_idx = find(strcmp(metadata.taskvars.sequence,'R'));
            metadata.taskvars.trial_iv_R.tstart = metadata.taskvars.trial_iv.tstart(right_idx);
            metadata.taskvars.trial_iv_R.tend = metadata.taskvars.trial_iv.tend(right_idx);
            
        end
        
%         if cfg.verbose
%             fprintf('start dt:\n'); disp(start_dt);
%             fprintf('end dt:\n'); disp(end_dt);
%         end
        
    case 'evt'
        
        error('Not yet implemented.');
        
        % load events
        
        % for each trial start, find next reward trigger, then preceding
        % photobeam; use those to restrict
        
    otherwise
        error('Unknown mode %s.',cfg.mode);
    
end





