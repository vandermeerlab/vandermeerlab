function [peth_out, all_trials] = tsdPETH_fast(cfg_in,tsd_in,t_in)
% function [peth_out, all_trials]  = tsdPETH(cfg_in,tsd_in,t_in)
%
% compute peri-event time histogram (average)
%
% INPUTS:
%
% tsd_in: tsd to compute average over
% t_in: times to average: can be raw timestamps, ts or iv data. If iv, cfg.window is ignored
%
% OUTPUTS:
%
% peth_out: tsd with PETH
% all_trials: tsd with individual trials, such that t_in(n) = all_trials.data(:, n)
%
% cfg options:
%
% cfg_def.window = [-2 2]; % start and end times of window (in s)
% cfg_def.dt = []; % time step, specify this for 'interp' mode
% cfg_def.mode = 'raw'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
% cfg_def.interp_mode = 'linear';
%
% EXAMPLE USAGE:
%
% cd('R042-2013-08-18'); LoadExpKeys; LoadMetadata;
% cfg = []; cfg.fc = ExpKeys.goodSWR(1);
% lfp = LoadCSC(cfg);
% peth = TSDpeth([],lfp,metadata.taskvars.trial_iv.tstart(2:end));
% plot(peth);
%
% NOTES:
%
% 'raw' mode collects tsd samples that fall within the specified intervals.
% Because it does not interpolate, it can easily fail if there are gaps in
% the data or the diffs between samples are unequal.
%
% 'interp' interpolates on a fixed timebase, so guarantees the number of
% samples will be constant.
%
% multidimensional TSDs are not yet supported!
%
% Manish 2024-04-07 (modified from TSDpeth.m)

cfg_def = [];
cfg_def.window = [-2 2]; % start and end times of window (in s)
cfg_def.dt = []; % time step, specify this for 'interp' mode
cfg_def.mode = 'raw'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
cfg_def.interp_mode = 'linear';

cfg = ProcessConfig2(cfg_def,cfg_in);

if ~CheckTSD(tsd_in)
    error('Incorrectly formed tsd input.');
end

if size(tsd_in.data,1) > 1
    error('More than 1 dimension in TSD input not yet supported.');
end

if ~isfield(t_in,'type') % assume raw times
    this_iv = iv(t_in + cfg.window(1), t_in + cfg.window(2));
else
    switch t_in.type
        case 'iv'
            
            if ~CheckIV(t_in) % this should check that there are no t_end before t_start, etc...
                error('Incorrectly formed iv input.');
            end
            
            this_iv = t_in;
               
        case 'ts'
            
            if length(t_in.t) ~= 1
                error('ts input must have exactly one .t cell');
            end
            
            this_iv = iv(t_in.t{1} + cfg.window(1), t_in.t{1} + cfg.window(2));
            
        otherwise
            
            error('t input must be raw times, iv or ts data')
        
    end
end

nT = length(this_iv.tstart);

switch cfg.mode
    
    case 'raw'
        
        % make big idx matrix
        start_idx = nearest_idx3(this_iv.tstart,tsd_in.tvec);
        end_idx = nearest_idx3(this_iv.tend,tsd_in.tvec);
        
        if length(unique(end_idx-start_idx)) ~= 1 % unequal length trials
           error('Raw mode requires equal tsd samples for each trial.');
        end

        % Number of elements in each sample
        num_el = end_idx(1) - start_idx(1) + 1; % guaranteed to be the same

        % start_idx needs to be a column matrix for this to work
        if isrow(start_idx)
            start_idx  = start_idx';
        end
        % Get a big matrix of indices to be extracted, extract all the
        % elements and then reshape
        big_data_idx = repmat(start_idx,1,num_el) + (0:num_el-1);
        big_data_idx = big_data_idx(:)';
        big_data = tsd_in.data(big_data_idx);
        out_data = reshape(big_data, [nT, num_el]);
        out_tvec = tsd_in.tvec(start_idx(1):end_idx(1));
        out_tvec = out_tvec-nanmean(out_tvec); % this is an approximation -- depends on exact spacing of input tsd
                
    case 'interp'
        
        if isempty(cfg.dt)
            error('interp mode requires cfg.dt to be specified.');
        end

        % First make a massive 1-D array with all the timestamps we will
        % require; 
        num_el = length(cfg.window(1):cfg.dt:cfg.window(2));
        
        % IMPORTANT: This assumes that this_iv.tstart is a column vector!
        big_tvec = repmat(this_iv.tstart,1,num_el) + cfg.dt*(0:num_el-1);
        big_tvec = big_tvec(:)'; % big_tvec must be row verctor
 
        % big_tvec is not necessarily in sorted order because the intervals
        % could be overlapping. So we first sort it, then run interp1,
        % transfer back to the unsorted version and then reshape!

        [~, sort_idx] = sort(big_tvec);
        [~, unsort_idx] = sort(sort_idx); % 

        big_data_sorted = interp1(tsd_in.tvec,tsd_in.data, big_tvec(sort_idx),cfg.interp_mode);
        % Unsort the data for reshaping
        big_data_og = big_data_sorted(unsort_idx);
        out_data = reshape(big_data_og, [nT,num_el]);
        out_tvec = cfg.window(1):cfg.dt:cfg.window(2);    
    otherwise
        
        error('unknown cfg.mode %s',cfg.mode);
end

% average and package
peth_out = tsd;
peth_out.data = nanmean(out_data);
peth_out.tvec = out_tvec;
CheckTSD(peth_out);

all_trials = tsd;
all_trials.data = out_data;
all_trials.tvec = out_tvec;
CheckTSD(all_trials);