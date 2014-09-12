function [MCC,redraw,rekey,undoable] = RemoveDoubleCounts(MCC,varargin)
% [MCC,redraw,rekey,undoable] = RemoveDoubleCounts(MCC,varargin)
%
% INPUTS
%    MCC
%
% OUTPUTS
%    updated MCC
%
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder

% MvdM 08
%
% remove all points which precede another point in the cluster by <1 ms
% and put them in a new cluster

redraw = true; rekey = true; undoable = true;

global MClust_FeatureTimestamps MClust_Clusters MClust_Hide

inCluster = FindInCluster(MCC);
timestamps = MClust_FeatureTimestamps(inCluster);
[ts_sorted,ts_ind] = sort(timestamps,'ascend'); % timestamps are in 10ths of ms, we think
ISIs = [Inf; diff(ts_sorted)];

badISIs = find(ISIs < 10); % less than 1 ms apart: double count
goodISIs = find(ISIs >= 10);

if ~isempty(badISIs)
    badSpikes = badISIs - 1; % delete first spike in doublets
    badSpikesInd = ts_ind(badSpikes);
    
    % create new cluster with bad spikes
    newCluster = mccluster('DoubleCounted');

    newCluster.myPoints = MCC.myPoints(badSpikesInd);
    newCluster.myOrigPoints = newCluster.myPoints;
    newCluster.name = 'DoubleCounted';
    
    MClust_Clusters{end+1} = newCluster;
	MClust_Hide(end+1) = 0;
    
    % delete bad spikes from cluster
	MCC.myOrigPoints(badSpikesInd) = [];
    MCC.myPoints(badSpikesInd) = [];   

end








