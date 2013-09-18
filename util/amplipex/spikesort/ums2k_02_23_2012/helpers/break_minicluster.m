function [spikes, new_clust] = break_minicluster(spikes, mini)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% break_minicluster - splits a mini cluster into two mini-clusters
%
% Usage:
%   [spikes, new_clust] = break_minicluster(spikes, mini)
%
% Description:  
%    Splits a mini-cluster into 2 mini-clusters by projecting the waveforms
% onto the mini-clusters's first principal component and dividing the 
% waveforms evenly between the 2 clusters. The function takes care of all
% book-keeping, such as resetting cluster labels and the cluster colors. 
%
% Input: 
%   spikes  - a spikes object
%   mini    - ID of the mini-cluster to split

% Output:
%   spikes    - the modified spikes object 
%   new_clust - the cluster ID for the new mini-cluster that was created 
% 

    %%%%%%%%%% ARGUMENT CHECKING
    if (~isfield(spikes.info, 'tree'))
        error('SS:hierarchical_information_unavailable', 'Hierarchical clustering must be performed before attempting to merge clusters.');
    end
    tree = spikes.info.tree;

    %  project waveforms in mini-cluster onto their principal components 
    members = find( spikes.info.kmeans.assigns == mini );
    proj = pcasvd(spikes.waveforms(members,:));

    % exactly half  are switched to the new cluster
    switching = members( find( proj(:,1) < median(proj(:,1)) ) );
    new_clust = spikes.info.kmeans.num_clusters + 1; % add 1 to the previously high ID
    spikes.info.kmeans.assigns( switching ) = new_clust;
    spikes.info.kmeans.num_clusters = new_clust;
    
    % add entry to beginning of aggregation tree that the new cluster was
    % added to the other
    tree = [mini new_clust; tree];
   
    % reset label on this cluster
    example = find( ismember( spikes.info.kmeans.assigns, mini ), 1); % sniff for ID by grabbing one of its members
    cluster = spikes.assigns( example );
    pos = find(spikes.labels(:,1) == cluster);
    spikes.labels( pos,2) = 1; 
    
    % and redo the color scheme
    spikes.info.kmeans.colors = emphasize_cluster_colors(spikes);