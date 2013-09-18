function [p,mu,stdev,n,x] = ss_undetected(spikes,use)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/201
%
% ss_undetected - wrapper for undetected.m
%
% Usage:
%     [p,mu,stdev,criteria] = ss_undetected(spikes,use)
%
% Description:  
%    Estimates how fraction of spikes associated with a given cluster that 
% may have gone undetected, by fitting a Gaussian with a missing tail to a 
% histogram of the detection criterion applied to each spike event.  
%
% See undetected.m for more information.
%
% Input:
%   spikes - spikes structure
%   use    - a cluster ID or an array describing which spikes to 
%            use as the first cluster
%          - see get_spike_indices.m
%
% Output:
%  p            - estimate of probability that a spike is missing because it didn't reach threshhold
%  mu           - mode estimated for distribution of minimum values
%  stdev        - standard deviation estimated for distribution of minimum values
%  n            - bin counts for histogram used to fit Gaussian
%  x            - bin centers for histogram used to fit Gaussian
%

    % apply detection criteria to each waveform on each channel in cluster
    select = get_spike_indices(spikes, use );      
    waveforms = spikes.waveforms(select,:,:);
    
    % threshes
    threshes = spikes.info.detect.thresh;

    criteria_func = spikes.params.detect_method;
    
    % call undetected
    [p,mu,stdev,n,x] = undetected(waveforms,threshes,criteria_func);
    
end







