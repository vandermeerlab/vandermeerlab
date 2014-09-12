function spikes =  reintegrate_outliers( spikes, indices, mini )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% reintegrate_outliers - handles logic of putting outlier spikes back
%                        into the main spikes object
%
% Usage:
%   spikes =  reintegrate_outliers( spikes, indices, subcluster )
%
% Description:
%    Handles logic of reintegrating outliers into the spikes object.  All
% information is copied back into the main spikes object and then deleted
% from the outlier subcluster.  The spikes are placed in a given 
% mini-cluster.
%
% Input: 
%   spikes - a spikes object
%   indices  - list of index into which outliers should be removed
%   mini    - ID of minicluster to add these spikes into
%           - either choose an existing minicluster or type 'new' to create
%             a new minicluster
%
% Output:
%   spikes - the modified spikes object
%

    % first check if there is an outliers object
    if ~isfield( spikes.info,'outliers'), error('No outliers present.');end
        
    if isequal( mini, 'new')
       mini = max( [spikes.info.outliers.subassigns spikes.info.kmeans.assigns] ) + 1; 
    end
    num_outliers = length(indices);

    % are we putting the clusters into a new or old minicluster?
    isnew = ~ismember( mini,spikes.info.kmeans.assigns);
    if ~isnew
        newassign =  spikes.assigns( find(spikes.info.kmeans.assigns,1) );
    else
        newassign =  mini; 
    end
    
    % add back into main spikes object
    o = spikes.info.outliers;
    spikes.waveforms = [spikes.waveforms; o.waveforms(indices,:,:)];
    spikes.info.kmeans.assigns = [spikes.info.kmeans.assigns mini*ones([1 num_outliers]) ];
    spikes.trials = [spikes.trials o.trials(indices)];
    spikes.spiketimes = [spikes.spiketimes o.spiketimes(indices)];
    spikes.unwrapped_times = [spikes.unwrapped_times o.unwrapped_times(indices)];
    spikes.assigns = [spikes.assigns newassign*ones([1 num_outliers]) ];
    spikes.info.pca.u = [spikes.info.pca.u; o.pca.u(indices,:) ];
    
    % remove from outliers
    spikes.info.outliers.waveforms(indices,:,:)=[];
    spikes.info.outliers.subassigns(indices)=[];
    spikes.info.outliers.spiketimes(indices) = [];
    spikes.info.outliers.unwrapped_times(indices) = [];
    spikes.info.outliers.trials(indices) = [];
    spikes.info.outliers.pca.u(indices,:) = [];
    
    % update labels
    % recolor if adding new cluster
    if isnew
        spikes.info.kmeans.colors = emphasize_cluster_colors( spikes );
      spikes.labels(end+1,:) = [newassign 1];
      spikes.info.kmeans.num_clusters = spikes.info.kmeans.num_clusters + 1; 
    else
      which = find(spikes.labels(:,1) == newassign );
      spikes.labels(which,:) = [newassign 1];
    end

end