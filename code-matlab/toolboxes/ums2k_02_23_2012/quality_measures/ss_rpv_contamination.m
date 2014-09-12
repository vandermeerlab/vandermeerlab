function [ev,lb,ub,RPV] = ss_rpv_contamination( spikes, use  )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% ss_rpv_contamination - wrapper for rpv_contamination.m
%
% Usage:
%   [lb,ub,ev] = ss_rpv_contamination(spikes,use)
%
% Description:
%   Estimates contamination of a cluster based on refractory period
% violations (RPVs).  See rp_violations.m for more information.  This 
% function simply allows rp_violations to be called simply on a SPIKES
% structure.
%
% Input: 
%   spikes - spikes structure
%   use    - a cluster ID or an array describing which spikes to 
%            use in analyzing refractory period violations
%          - see get_spike_indices.m
%
% Output:
%   ev   - expected value of % contamination,
%   lb   - lower bound on % contamination, using 95% confidence interval
%   ub   - upper bound on % contamination, using 95% confidence interval
%   RPV  - number of refractory period violations observed in this cluster
%

    % get spiketime data
    select = get_spike_indices(spikes, use );      
    spiketimes = sort( spikes.unwrapped_times(select) );

    % get parameters for calling rp_violations
    N = length(select);
    T = sum( spikes.info.detect.dur );
    RP = (spikes.params.refractory_period - spikes.params.shadow) * .001; 
    RPV  = sum( diff(spiketimes)  <= (spikes.params.refractory_period * .001) );

    % calculate contamination
    [ev,lb,ub] = rpv_contamination(N, T, RP, RPV );
  
end
