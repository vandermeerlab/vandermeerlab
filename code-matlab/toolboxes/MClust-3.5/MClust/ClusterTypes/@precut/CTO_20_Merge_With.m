function [MCC, redraw, rekey, undoable] = CTO_MergeWith(MCC, varargin)

% f = AddPoints(MCC, varargin)
%
% INPUTS
%     MCC - a MCCluster
% needs iClust from parent
%       figHandle
%
% OUTPUTS
%     MCC - The updated cluster
%
% 
% ncst 26 Nov 02
% ADR 2008
% precut
%

redraw= true; rekey = true; undoable = true;

global MClust_Clusters MClust_Colors
global MClust_Hide

extract_varargin;

mergeSet = inputdlg(['Merge Cluster ' num2str(iClust) ' with ...? (Enter a single cluster or a blank separated list of clusters...)']);
if isempty(mergeSet), return, end
mergeSet = str2num(mergeSet{1}); %#ok<ST2NM>

newCluster = precut('MergedCluster');

newCluster.myPoints = FindInCluster(MCC);
newCluster.name = ['Merge of (' num2str(iClust)];

mergeFail = 0;
for iMC = 1:length(mergeSet)
	if ~mergeFail && mergeSet(iMC) > 0 && mergeSet(iMC) <= length(MClust_Clusters) && mergeSet(iMC) ~= iClust
		newCluster.myPoints = cat(1, newCluster.myPoints, FindInCluster(MClust_Clusters{mergeSet(iMC)}));
		newCluster.name = cat(2, newCluster.name, ',', num2str(mergeSet(iMC)));
	else
		mergeFail = mergeSet(iMC);
	end
end

if ~mergeFail
	newCluster.name = cat(2, newCluster.name, ')');
	newCluster.myPoints = unique(newCluster.myPoints);
	
	MClust_Clusters{end+1} = newCluster;
	MClust_Hide(end+1) = 0;
    MClust_Colors(length(MClust_Hide),:) = MClust_Colors(iClust+1,:)*0.95;

	warndlg(['Cluster ' num2str(iClust) ' merged with cluster(s) ' ...
		num2str(mergeSet) ' into new cluster ' num2str(length(MClust_Clusters)) '.'], 'Merge successful');
else
	errordlg(['Cannot merge clusters ' num2str(iClust) ' and ' num2str(mergeFail) '.']);
end
