function RunKKwik_Callbacks

% CBP_Callbacks
%
% Add a limit on the best projection of two clusters.
% 
% If cluster 0 is selected as the second cluster, then the program attempts
% to separate the cluster from the noise.  Takes a long time if you have a
% lot of spikes.
%
% INPUTS
% 
% NONE
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust folder
% Also requires Find_Best_Projection.m

% NCST 2003
%

    cboHandle = gcbo;
    figHandle = gcf;
    callbackTag = get(cboHandle, 'Tag');
  
global MClust_FeatureNames MClust_FeatureSources MClust_FDfn 
global MClust_Clusters MClust_Colors

% get parameters from figure
newFeatures = get(findobj(figHandle,'Tag','FeaturesUseListbox','TooltipString', ...
    'These features will be used for cluster separation.'),'String');

iClust = get(findobj('Tag','RunKlustaKwik'), 'UserData');
minClusters = str2double(get(findobj(figHandle,'Tag','RunKlustaKwikMinClust'),'String'));
maxClusters = str2double(get(findobj(figHandle,'Tag','RunKlustaKwikMaxClust'),'String'));
otherParms = get(findobj(figHandle, 'Tag', 'OtherParameters'), 'String');

maxPossibleClusters = maxClusters;

close(figHandle);

% find KlustaKwik
KlustaKwikPath = which('KlustaKwik.exe');
if ~isempty(KlustaKwikPath)
    disp(['Using ' KlustaKwikPath]);
else
    disp('Did not find KlustaKwik.exe.');
    return
end     


% get features, construct FeatureData
[spikeIndex MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});

if isempty(newFeatures)
    errordlg('No features chosen', 'RunKKwik Error', 'modal');
    return;
else
	nFeatures = length(newFeatures);
	FeatureData = nan(length(spikeIndex), nFeatures);
	for iF = 1:nFeatures
		FeatureID = strmatch(newFeatures{iF}, MClust_FeatureNames);
		if isempty(FeatureID)
			errordlg('No matching feature found', 'RunKKwik Error', 'modal');
			return;
		elseif length(FeatureID) > 1
			errordlg('Too many matching features found', 'RunKKwik Error', 'modal');
			return;
		else
			temp = load(MClust_FeatureSources{FeatureID,1}, '-mat', 'FeatureData');
			FeatureData(:,iF) = temp.FeatureData(spikeIndex,MClust_FeatureSources{FeatureID,2});
		end
	end
end

% Set up parameters
[fpath fname fext] = fileparts(MClust_FDfn);
FETfn = fullfile(fpath,[fname '-Cluster' num2str(iClust)]);
WriteFeatureData2TextFile(FETfn,FeatureData);

% GO
file_no = 1;
parameter_string = ['-MinClusters ' num2str(minClusters) ...
        ' -MaxClusters ' num2str(maxClusters) ...
        ' -MaxPossibleClusters ' num2str(maxPossibleClusters) ...
		' ', otherParms, ' ', ...
		' -UseFeatures ' num2str(repmat(1,size(FeatureData,2),1))'];
    % added " " to handle Windows spaces
COMMAND = ['! "' KlustaKwikPath '" "' FETfn '" ' num2str(file_no) ' ' parameter_string ];
disp(FETfn);
disp(['Number of spikes: ' num2str(size(FeatureData,1))]);
%                disp(['Estimated time required to run KlustaKwik.exe on this file is ' num2str(EstDuration(i)*60) ' seconds (or ' num2str(EstDuration(i)/60) 'hours)']);
COMD_output = evalc(COMMAND);
% Find the output that you are going to display
COMD_filebase = findstr(COMD_output,'FileBase');
COMD_dim = findstr(COMD_output,'dimension');
COMD_time = findstr(COMD_output,'That took');

try
    disp(COMD_output([COMD_filebase:COMD_dim + 13,COMD_time:end]));
catch
    disp(['Did not find output that was searched for, all output shown:  ' COMD_output]);
end

% diary on
% disp(COMD_output)
% diary off

disp(' ')

clusters = load([FETfn '.clu.1']);
clusters = clusters(2:end);            
unique_clusters = unique(clusters);
n_clusters = length(unique_clusters);

cmap = colormap('lines');
for iC = 1:n_clusters
    idx = logical(clusters == unique_clusters(iC));
	nC = length(MClust_Clusters);
	MClust_Clusters{nC + 1} = precut(['Cluster ' num2str(iClust) '_sub' num2str(iC)]);
	MClust_Clusters{nC + 1} = AddIndices(MClust_Clusters{nC+1}, spikeIndex(idx));
	if size(MClust_Colors,1) < nC+2
		MClust_Colors(nC + 2,:) = (MClust_Colors(iClust,:) + cmap(iC,:))/2;
	end
end

figHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');
MClustCutterClearClusterKeys(figHandle);
MClustCutterRedrawClusterKeys(figHandle, max(0,length(MClust_Clusters)-16));
MClustCutterCallbacks('RedrawAxes')

% ClusterCutWindow = findobj('Type','figure','Tag', 'ClusterCutWindow');
% if ~isempty(ClusterCutWindow)
%     close(ClusterCutWindow);
% end
% 
% CHDrawingAxisWindow = findobj('Type','figure','Tag', 'CHDrawingAxisWindow');
% if ~isempty(CHDrawingAxisWindow)
%     close(CHDrawingAxisWindow);
% end

fn = FindFiles([FETfn '*'], 'CheckSubdirs', 0); % ADR 8 Aug 2008
for iFN = 1:length(fn)
    eval(['! del ' fn{iFN}]);
end


%===============================================================================
function WriteFeatureData2TextFile(file_name, FeatureData)
%
% write featuredata from memory to a text file for input into KlustaKwick.exe
%
file_no = 1;
fid = fopen([ file_name '.fet.' num2str(file_no)],'w');
[n_points, n_features] = size(FeatureData);
fprintf(fid,'%3d \n',n_features);
for ii = 1:n_points
    fprintf(fid,'%f\t',FeatureData(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
