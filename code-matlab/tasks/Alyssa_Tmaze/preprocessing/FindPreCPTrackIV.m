function PreCPTrackIV = FindPreCPTrackIV(cfg_in,pos,metadata,ExpKeys)
% function PreCPTrackIV = FindPreCPTrackIV(cfg_in,pos,metadata,ExpKeys)
%
% Finds MotivationalT intervals where animal is on the track but before
% the choice point
%
% INPUTS:
%
% *pos: position data in pixels (not cm!), obtained from LoadPos([])
% *metadata, ExpKeys
%
% OUTPUTS:
%
% *PreCPTrackIV: intervals where animal is on the track but before
% the choice point
%
% CONFIG DEFAULTS:
% 
% cfg_def.startbuffer = Inf; % padding (s) from end of preceding pedestal epoch to allow for transit time
% NOTE, output will use metadata.taskvars.trial_iv.tstart if padding time ends up later than metadata.taskvars.trial_iv.tstart
% cfg_def.endbuffer = 0; % padding (s) before end of run before reaching CP (note: -1 means stop trial 1 s before reaching CP!)
% cfg_def.verbose = 1; % show some debugging output
%
% HOW HAS THIS BEEN TESTED?
%
% I looked at debugging output for all sessions (MvdM), no outliers found
% I plotted the position data restricted by the output for an example session, looked correct (MvdM)
%
% ANY CAVEATS?
%
% This works by finding the sample withe smallest distance to the CP, and restricts to all times before that.
% If the animal passes the CP multiple times, there is no guarantee that this approach will pick up the first pass.
%
% MvdM 2017

cfg_def.startbuffer = Inf; % padding (s) from end of preceding pedestal epoch to allow for transit time
% NOTE, output will use metadata.taskvars.trial_iv.tstart if padding time ends up later than metadata.taskvars.trial_iv.tstart
cfg_def.endbuffer = 0; % padding (s) before end of run before reaching CP (note: -1 means stop trial 1 s before reaching CP!)
cfg_def.verbose = 1; % show some debugging output
cfg = ProcessConfig(cfg_def,cfg_in);

% initialize empty iv
PreCPTrackIV = iv;

% get all times identified as end of a pedestal epoch; this includes the prerecord
start_values = cat(1,ExpKeys.prerecord(2),metadata.taskvars.rest_iv.tend);

for iT = 1:length(metadata.taskvars.sequence) % loop over trials
      
    % find pedestal end time immediately preceding this trial start
    idx = nearest_idx3(metadata.taskvars.trial_iv.tstart(iT),start_values,-1); % -1: look backwards
    this_start = start_values(idx);
    
    if this_start + cfg.startbuffer > metadata.taskvars.trial_iv.tstart(iT)
        % if somehow the trial starts before the end of the preceding pedestal epoch plus buffer,
        % then the time identified as trial start in metadata takes precedence. This can happen
        % if the experimenter and rat were faster to start the next trial than the time in 
        % cfg.buffer.
        this_start = metadata.taskvars.trial_iv.tstart(iT);
    end
    
    this_end = metadata.taskvars.trial_iv.tend(iT); % this is the original end time, at the end of maze arms
    
    %
    this_cp = metadata.coord.chp;
    
    this_pos = restrict(pos,this_start,this_end); % position data for this trial only
    this_pos.data(3,:) = sqrt((this_pos.data(1,:) - this_cp(1)).^2 + (this_pos.data(2,:) - this_cp(2)).^2); % distance to cp
    
    [min_dist,min_dist_idx] = min(this_pos.data(3,:));
    
    if cfg.verbose
       fprintf('FindPreCPTrackIV.m: Trial %d, minimum distance from CP is %.1f pix\n',iT,min_dist);
    end
    
    this_end = this_pos.tvec(min_dist_idx) + cfg.endbuffer;
    
    PreCPTrackIV.tstart(iT) = this_start;
    PreCPTrackIV.tend(iT) = this_end;
    
end