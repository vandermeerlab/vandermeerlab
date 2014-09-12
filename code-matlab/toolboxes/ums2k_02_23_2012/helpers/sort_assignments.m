function reassignments = sort_assignments(assignments)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/08/2010
%
% sort_assignments - renumber cluster identities based on size of cluster
%
% Usage:
%   reassignments = sort_assignments(assignments)
%
% Description:
%   Takes a list of assignment numbers and reassigns label numbers such
% that the largest size group is assigned label '1', the next largest is 
% assigned label '2', and so on.
%
% Input: 
%   assignments - [1 X N] list of cluster ID for each spike event
%
% Output:
%   reassignments - [1 X N] list of cluster ID with new IDs
%

clusters = unique(assignments);  % get a list of unique labels . . .
numclusts = length(clusters);    % 
clustsize = zeros(numclusts,1);  %
for clust = 1:numclusts          % ... and count # elements assigned to each label
    clustsize(clust) = length(find(assignments == clusters(clust)));
end

% create a matrix with cols [old_label  num_elements] and (descending) sort on num_elemebts
reassign_list = flipud(sortrows([clusters, clustsize], 2));

% . . . and use that table to translate the original assignment list
reassignments = zeros(size(assignments));
for clust = 1:numclusts
    reassignments(find(assignments == reassign_list(clust,1))) = clust;
end

