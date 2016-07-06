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
%
% MvdM 2016-05-17 initial version

cfg_def = [];
cfg_def.maxJump = 20; % number of bins decoded Z can jump without breaking sequence
cfg_def.minLength = 3;

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
prev_val = -100;

tstart = []; tend = [];
for iT = 1:length(decoded_z.data)
    
    this_jump = abs(decoded_z.data(iT)-prev_val);
    
    if this_jump < cfg.maxJump
        
        % check if sequence started
        if isempty(this_start_idx)
           this_start_idx = iT-1; 
        end
        
    else % sequence broken
    
        % check if long enough
        if iT - this_start_idx >= cfg.minLength % criterion met
            nDetected = nDetected + 1;
            tstart(nDetected) = decoded_z.tvec(this_start_idx);
            tend(nDetected) = decoded_z.tvec(iT-1);
        end
        this_start_idx = [];
        
    end
    
    prev_val = decoded_z.data(iT);
    
end

seq_iv = iv(tstart,tend);