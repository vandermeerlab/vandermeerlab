function confusion = ss_gaussian_overlap( spikes, use1, use2 )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_gaussian_overlap - wrapper for gaussian_overlap.m
%
% Usage:
%     C = ss_gaussian_overlap( spikes, use1, use2 )
%
% Description:  
%   Estimates the overlap between 2 spike clusters by fitting with two
% multivariate Gaussians.  See gaussian_overlap.m for more information.  This 
% function simply allows gaussian_overlap to be called simply on a SPIKES
% structure.
%
% Input: 
%   spikes - spikes structure
%   use1    - a cluster ID or an array describing which spikes to 
%            use as the first cluster
%   use2    - a cluster ID or an array describing which spikes to
%            use as the second cluster
%           - see get_spike_indices.m
%
% Output:
%   C   - a confusion matrix
%   C(1,1) - False positive fraction in cluster 1 (waveforms of neuron 2 that were assigned to neuron 1)
%   C(1,2) - False negative fraction in cluster 1 (waveforms of neuron 1 that were assigned to neuron 2)
%   C(2,1) - False negative fraction in cluster 2 
%   C(2,2) - False positive fraction in cluster 2
%
 
    % get both sets of waveforms
    select1 = get_spike_indices(spikes, use1 );      
    waveforms1 = spikes.waveforms( select1, : );
    select2 = get_spike_indices(spikes, use2 );      
    waveforms2 = spikes.waveforms( select2, : );
   
    % call gaussian_overlap
    confusion = gaussian_overlap( waveforms1, waveforms2 );

end

