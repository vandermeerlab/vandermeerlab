function colors = emphasize_cluster_colors(spikes) 
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% emphasize_cluster_colors - make color matrix for each minicluster ID 
%                            where the colors for clusters are separated
%
% Usage:
%      colors = emphasize_cluster_colors(spikes) 
%
% Description:  
%    Creates a color matrix with one row of RGB values for each
% minicluster.  The entire color spectrum is used, but the entries for 
% major clusters are assigned first so that they can be well separated
% in color space.
%
% Input: 
%   spikes - a spikes object
%
% Output:
%   colors - [3 x N] color map where N = no. of mini-clusters
%

    % create a large enough color map
    clust_list = sort(unique(spikes.assigns));
    num_subclusters = spikes.info.kmeans.num_clusters;
    cmap       = jetm(  num_subclusters );
    
    % choose cluster indices to be well spread out
    cluster_indices = round( linspace(1,num_subclusters,length(clust_list) ) );
    colors     = zeros(size(cmap) );
    colors(clust_list,:) = cmap(cluster_indices,:);
    
    % randomly assign the rest of the colors
    remaining_indices = setdiff(1:num_subclusters, cluster_indices );
    remaining_list    = setdiff(1:num_subclusters, clust_list );
    colors( remaining_list,: ) = cmap( randperm(length(remaining_indices)), : );
