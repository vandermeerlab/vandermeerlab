function c = ss_censored( spikes, use )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_censored - wrapper for censored.m
%
% Usage:
%           c = ss_censored( spikes, use )
%
% Description:  
%   Estimates the percent of a cluster that is lost due to censoring.  
% See censored.m for more information.  This function simply allows 
% censored to be called more directly on a SPIKES structure.
%
% Input: 
%   spikes - spikes structure
%   use    - a cluster ID or an array describing which spikes to 
%            use in analyzing censoring
%          - see get_spike_indices.m
%
% Output:
%       c    - estimated false negative fraction from censoring
%
  
    % determine number of events outside cluster of interest
    select = get_spike_indices(spikes, use );
    M      = length(spikes.spiketimes) - length(select);
    if isfield(spikes.info,'outliers') % if there are outliers, include them
        M = M + length(spikes.info.outliers.spiketimes);
    end
    
    % get other parameters
     tau_c = spikes.params.shadow;
     T     =  sum( spikes.info.detect.dur );

    % finally, calculate censored period
    c = censored( tau_c, M, T );
    