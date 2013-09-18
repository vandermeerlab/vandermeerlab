function spikes = merge_clusters(spikes, to, from)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% merge_clusters - Merge two clusters after automatic hierarchical aggregation.
%
% Usage:
%       spikes = merge_clusters(spikes, to, from)
%
% Description:
%    SPIKES = MERGE_CLUSTERS(SPIKES, TO, FROM) takes and returns a spike-
%    sorting object SPIKES.  SPIKES must have gone through a hierarchical
%    clustering aggregation (e.g., SS_AGGREGATE) previous to this function call.
%   
%    All spikes belonging to the cluster whose label number is given by FROM
%    are merged into the cluster with label TO.  The hierarchical clustering
%    assignments and aggregation tree are modified to reflect this change.
%    If SPIKES contains a SPIKETIMES vector, the ratio of the interspike  
%    interval count below 2 msec to the count below 10 msec is computed and
%    entered into the SPIKES.HIERARCHY.TREE.
%

%%%%%%%%%% ARGUMENT CHECKING
if (~isfield(spikes.info, 'tree'))
	error('SS:hierarchical_information_unavailable', 'Hierarchical clustering must be performed before attempting to merge clusters.');
elseif (~all(ismember([to,from], unique(spikes.assigns))))
    error('SS:cluster_numbers_not_found', 'One or both of the cluster labels supplied does not exist.');
elseif ((length(from) > 1) | (length(to) > 1))
    error('SS:one_at_a_time', 'The ''from'' and ''to'' labels must be scalar values.');
end

%%%%%%%%%% MERGING ASSIGNMENTS
members_from = find(spikes.assigns == from);    % Get list of spikes to move ...
orig_members_to = find(spikes.assigns == to);   %   (we need this below)
spikes.assigns(members_from) = to;              % ... and relabel them.
spikes.info.tree = [spikes.info.tree; to from];  % add entry to tree structure

spikes.info.kmeans.colors = emphasize_cluster_colors(spikes); % recalculate colors 

% reset labels to default value (1)
spikes.labels( spikes.labels(:,1) == from,: )  = [];
spikes.labels( spikes.labels(:,1) == to, 2 )  = 1;
