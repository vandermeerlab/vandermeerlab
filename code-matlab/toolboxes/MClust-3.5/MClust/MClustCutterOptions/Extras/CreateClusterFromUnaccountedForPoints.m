function [redraw, rekey, undoable] = CreateClusterFromUnaccountedForPoints

% [redraw, rekey, undoable] = CreateClusterFromUnaccountedForPoints
%
% INPUTS
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
%
% ADR 2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

redraw = true; rekey = true; undoable = true;

global MClust_Clusters MClust_Hide MClust_Colors MClust_CurrentFeatureData
global MClust_AvailableClusterTypes

if isempty(MClust_Clusters)
	warndlg('No clusters to count.');
	return
end
MCC = precut('Unaccounted For Points');
MClust_ClusterIndex = ProcessClusters(MClust_CurrentFeatureData, MClust_Clusters);
f = find(MClust_ClusterIndex == 0);
MCC = AddIndices(MCC, f);

clusterType = MClust_AvailableClusterTypes{get(findobj('Tag', 'AddAsType'), 'value')};
MCC = feval(clusterType, 'Unaccounted For Points', MCC);

MClust_Clusters{end+1} = MCC;
MClust_Hide(length(MClust_Clusters)+1) = 0;
MClust_Colors(length(MClust_Clusters)+1,:) = [0.75 0.75 0.75];
  
