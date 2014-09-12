function spikes = remove_outliers( spikes, which )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% remove_outliers - handles logic of labeling set of spikes as outliers
%
% Usage:
%   spikes = remove_outliers(spikes, bad  )
%
% Description:
%    Handles logic of labeling set of spikes as outliers.  Instead of 
% throwing away outliers, we place them in a sub-structure:
% spikes.info.outliers.  Here we store all information associated with
% the outlier spike for later examination or re-integration.  Any cluster
% that has had outliers removed has its label reset.
%
% Input: 
%   spikes - a spikes object
%   which  - specifies the spikes to be removed -- see get_spike_indices.m 
%
% Output:
%   spikes - the modified spikes object
%

    % check inputs
    if ~isfield(spikes,'assigns'), error( 'Outlier removal must be performed after cluster aggregation.'); end
    
    % make lists
    subcluster_list = unique(spikes.info.kmeans.assigns);
    cluster_list    = unique(spikes.assigns);
    
    bad = get_spike_indices(spikes, which );

    if ~isfield(spikes.info,'outliers')
        
        % if we are dealing with outliers for the first time, the outliers
        % sub-structure needs to be initialized
        spikes.info.outliers.waveforms = [];
        spikes.info.outliers.subassigns = [];
        spikes.info.outliers.trials = [];
        spikes.info.outliers.spiketimes = [];
        spikes.info.outliers.pca.u = [];
        spikes.info.outliers.unwrapped_times = [];
        
    end
       
       % append data to outliers structure for associated fields
       spikes.info.outliers.waveforms = [spikes.info.outliers.waveforms; spikes.waveforms(bad,:,:) ];
       spikes.info.outliers.trials = [spikes.info.outliers.trials spikes.trials(bad) ];
       spikes.info.outliers.spiketimes = [spikes.info.outliers.spiketimes spikes.spiketimes(bad) ];
       spikes.info.outliers.subassigns = [spikes.info.outliers.subassigns spikes.info.kmeans.assigns(bad) ];
       spikes.info.outliers.pca.u = [spikes.info.outliers.pca.u; spikes.info.pca.u(bad,:) ];
      spikes.info.outliers.unwrapped_times = [spikes.info.outliers.unwrapped_times spikes.unwrapped_times(bad) ];

       % update tree if whole subclusters have been removed
       tempassigns = spikes.info.kmeans.assigns;  tempassigns(bad) = [];
       subcluster_list = setdiff(subcluster_list, unique(tempassigns) );
       for j = subcluster_list(1:end-1)
           spikes = split_cluster(spikes, j);           
       end

       % remove data from main spikes object
       spikes.waveforms(bad,:,:) = [];
       spikes.trials(bad) = [];
       spikes.spiketimes(bad) = [];
       spikes.unwrapped_times(bad) = [];
       spikes.info.kmeans.assigns(bad) = [];
       spikes.assigns(bad)=[];
       spikes.info.pca.u(bad,:) = [];

       % update label
       cluster_list = setdiff(cluster_list, unique(spikes.assigns) );
       which = ismember( spikes.labels(:,1), cluster_list );
       spikes.labels(which,:) = [];
    