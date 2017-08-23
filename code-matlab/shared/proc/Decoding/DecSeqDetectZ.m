function [seq_iv,decoded_z] = DecSeqDetectZ(cfg_in,P)
% function [seq_iv,decoded_z] = DecSeqDetectZ(cfg_in,P)
%
% detects intervals containing sequences of decoded locations
%
% INPUTS:
%
% P: tsd containing decoded posterior - output of DecodeZ()
%
% OUTPUTS:
%
% seq_iv: iv with start and end times of detected sequences
% decoded_z: tsd with locations of mode of decoded posterior (bin with
%  largest probability)
%
% CONFIGS:
%
% cfg_def.maxJump = 20; % number of space bins decoded Z can jump without breaking sequence
% cfg_def.minLength = 3; % number of time bins without break for sequence
%  to be included
% cfg_def.nMaxNanSkipSequential = 1; % number of consecutive NaN bins that can be skipped 
%   before sequence is terminated
% cfg_def.nMaxNanSkipTotal = Inf; % if non-integer, maximum fraction of sequence length 
%  that is NaN acceptable for sequence to be included; if integer, maximum NaN count
%
% MvdM 2016-05-17 initial version
% MvdM 2017-08-21 update to allow skips

cfg_def = [];
cfg_def.maxJump = 20; % number of bins decoded Z can jump without breaking sequence
cfg_def.minLength = 3; % number of time bins without break for sequence
cfg_def.nMaxNanSkipSequential = 0; % number of consecutive NaN bins that can be skipped 
% before sequence is terminated
cfg_def.nMaxNanSkipTotal = Inf; % if non-integer, maximum fraction of sequence length 
% that is NaN acceptable for sequence to be included; if integer, maximum NaN count

cfg = ProcessConfig(cfg_def,cfg_in);

% extract Zbin with largest posterior probability ("decoded location")
[~,decoded_z] = max(P.data);
decoded_z = tsd(P.tvec,decoded_z);

% set decoded z for nanbins to nan
nanbin = all(isnan(P.data));
decoded_z.data(nanbin) = nan;

% detect
nT = length(decoded_z.data);
nDetected = 0;

this_start_idx = [];
prev_val = -Inf;

tstart = []; tend = [];
current_skip = 0; total_skip = 0;
for iT = 1:length(decoded_z.data)
    
    sequence_broken = 0;
       
    if isnan(decoded_z.data(iT))
    
        % NaNs cannot start a sequence, but if sequence already started, may tolerate
        if ~isempty(this_start_idx)
            current_skip = current_skip + 1;
            total_skip = total_skip + 1;
            if current_skip > cfg.nMaxNanSkipSequential
                sequence_broken = 1; prev_val = NaN;
            end
        else
            prev_val = NaN;
        end
           
    end % could be else?
    
    this_jump = abs(decoded_z.data(iT)-prev_val);
        
    if this_jump <= cfg.maxJump % good sample 
        
        % check if sequence started; if not, start it
        if isempty(this_start_idx)
           this_start_idx = iT-1;
           total_skip = 0;
        end
        
        % also reset skip counter
        current_skip = 0;
    
    end
    
    if ~isnan(decoded_z.data(iT)) & this_jump > cfg.maxJump
       sequence_broken = 1; 
    end
        
    end_of_data = (iT == length(decoded_z.data));
    if sequence_broken | end_of_data % sequence broken, or end of data
    
        % check if long enough
        this_length = iT - this_start_idx - current_skip; % subtract trailing NaNs from length
        
        if end_of_data & ~sequence_broken
           this_length = this_length + 1; 
        end
        
        if this_length >= cfg.minLength % length criterion met
            sequence_pass = 1;
            
            % test for total NaNs here
            if mod(cfg.nMaxNanSkipTotal,1) == 0 % integer, parameter is a count
                if total_skip > cfg.nMaxNanSkipTotal
                   sequence_pass = 0; 
                end
            else % non-integer, parameter is a fraction of length
                if total_skip > floor(this_length.*cfg.nMaxNanSkipTotal)
                    sequence_pass = 0;
                end
            end
            
            if sequence_pass
                % if pass, store for output
                nDetected = nDetected + 1;
                tstart(nDetected) = decoded_z.tvec(this_start_idx);
                
                %if current_skip > 0
                %    nBack = current_skip - 1; % sequence can't end in NaN
                %else
                %    nBack = 0;
                %end
                if isnan(decoded_z.data(iT)), current_skip = current_skip - 1; end 
                
                % annoying edge case if last sample was part of a sequence
                if (iT == length(decoded_z.data)) & ~sequence_broken
                    lastSampleAdjustment = 1;
                else
                    lastSampleAdjustment = 0;
                end
                
                %tend(nDetected) = decoded_z.tvec(iT-1-nBack+lastSampleAdjustment); % fails with 1 trailing NaN --> isn't removed
                tend(nDetected) = decoded_z.tvec(iT-1-current_skip+lastSampleAdjustment); % fails with 2 trailing NaNs --> cuts off last correct point
            end
        end % of sufficient length
        this_start_idx = [];
      
    end % of sequence_broken
    
    if ~isnan(decoded_z.data(iT))
        prev_val = decoded_z.data(iT);
    end
    
end

seq_iv = iv(tstart,tend);