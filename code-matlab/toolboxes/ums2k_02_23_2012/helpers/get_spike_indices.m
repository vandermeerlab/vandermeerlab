function which = get_spike_indices(spikes, show )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% get_spike_indices - get the indices of the specified spikes
%
% Usage:
%   which = get_spike_indices(spikes, show )
%
% Description:
%   Gets the indices of a set of spikes in the spikes object.  The input
% show is flexible to allow spike selection in a variety of ways.
%
% Input: 
%   spikes - a spikes object
%   show   - Can be:   'all' - all spikes will be selected
%                       [1 x N] binary array - selects spikes where show contains the value 1
%                       [1 x M] index array where M < N - returns spikes that are members of 
%                                                         clusters in show subclusters if 
%                                                         clusters is unavailable)
%
% Output:
%   which  - [1 x K] list of indices for events specified in show 
%

    num_spikes = size(spikes.waveforms,1);

    % We turn show into a binary array
    
    if isequal(show,'all') % if selecting all, all values of show are 1
        show = ones( [1 num_spikes ] );

    elseif length(show) < num_spikes % if selecting by cluster ID
        
            if isfield( spikes,'assigns') % use cluster ID if available
                show = ismember( spikes.assigns, show );
                
            elseif isfield(spikes.info,'kmeans') % else use the subcluster IDs
                show = ismember( spikes.info.kmeans.assigns, show );
                
            else                % else returns all spikes
                show = ones( [1 num_spikes ] );
            end
    end

    which = find(show); % return events specified by show

