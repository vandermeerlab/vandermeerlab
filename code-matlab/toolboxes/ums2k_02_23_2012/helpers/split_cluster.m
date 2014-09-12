function [spikes, new_clust] = split_cluster(spikes, subclusters)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% split_cluster - split a cluster into two parts
%
% Usage:
%    [spikes, new_clust] = split_cluster(spikes, target)
%
% Description:  
%   Breaks off the specified sub-tree from the rest of the 
% aggregation hierarchy.  The sub-tree is specified by the list of all
% the mini-clusters that it contains.  The function takes care of all
% book-keeping, such as resetting cluster labels and the cluster colors. 
% NOTE: It is the calling function's responsibility to make sure that the 
% list of mini-clusters actually constitutes a valid sub-tree. 
%
% Input: 
%   spikes  - a spikes object
%   target  - an array of mini-clusters that constitute this sub-tree
%
% Output:
%   spikes    - the modified spikes object 
%   new_clust - the cluster ID for the new cluster that was created 
%
 
    %%%%%%%%%% ARGUMENT CHECKING
    if (~isfield(spikes.info, 'tree'))
        error('SS:hierarchical_information_unavailable', 'Hierarchical clustering must be performed before attempting to merge clusters.');
    end
    tree = spikes.info.tree;
    
    % Get ID for overall cluster  
    example = find( ismember( spikes.info.kmeans.assigns, subclusters(1) ), 1); % sniff for ID by grabbing one of its members
    cluster = spikes.assigns( example );

    % find any lines in the tree matrix involving these subclusters
    to   = ismember( tree(:,1), subclusters ); % operations where something was merged TO this sub-tree
    from = ismember( tree(:,2), subclusters ); % operation where something FROM this subtree was merged to somethign else
    
    % There are 2 situations we need to deal with.
    % SITUATION #1:  Nothing was added to this sub-cluster after it was fully formed.
    if ~any( to & ~from)
        tree( find( ~to & from ), : ) = []; % just deleted the entry for merging the subtree into something else
    
    % SITUATION #2: Other events were added to the sub-cluster after it was fully formed
    else
        % find the first time something else was added to our subcluster
        first = find( to & ~from,1 );
        alt_id    = tree(first,2);
        my_id     = tree(first,1);
 
        % delete this line
        tree(first,:) = [];
        
        % replace any mention of the tree
        idx = [first:size(tree,1)];
        i  = find(  tree(idx,1) == my_id );
        tree(first+i-1,1) = alt_id;
        i  = find(  tree(idx,2) == my_id );
        tree(first+i-1,2) = alt_id;
        
    end
 
  %
  % BOOK-KEEPING
  %
    
   % Now the tree is fine, but we need to make sure the IDs are correct, so
   % we re-number all events based on the merge operations in the tree.
   spikes.info.tree =tree;    
   spikes.assigns = spikes.info.kmeans.assigns;
   for j = 1:size(tree,1)
    spikes.assigns(  spikes.assigns == tree(j,2)) = tree(j,1);
   end
   
   % determine the new tree's ID -- it's the one that doesn't have a label yet
    new_clust = setdiff( unique(spikes.assigns), spikes.labels(:,1)' );
  
    % make the default label for the new cluster
    pos = find( spikes.labels(:,1) == new_clust); 
    if isempty(pos), pos = size(spikes.labels,1) + 1; end
      spikes.labels(pos,:) = [new_clust 1];
      
    % reset label for original cluster as well
    orig = find(spikes.labels(:,1) == cluster );
    spikes.labels( orig,2) = 1; 
    [junk,a] = sort(spikes.labels(:,1)); spikes.labels = spikes.labels(a,:);
    
    % recolor all clusters
    spikes.info.kmeans.colors = emphasize_cluster_colors(spikes);
    end    