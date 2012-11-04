function ViewAllDimensionPairs()

% ViewAllDimensionPairs

% find all feature names
global MClust_FeatureNames
nFeatures = length(MClust_FeatureNames); % all but time

ppL = nFeatures-1;
nL = nFeatures-1;

figure('NumberTitle', 'off', 'Name', 'FeaturePairs');

iPlot = 0;
for iC = 2:nFeatures
	for jC = 1:(nFeatures-1)
		iPlot = iPlot+1;
		if jC < iC
			subplot(ppL, nL, iPlot);
			hold on
			DrawAxes(iC, jC);
			xlabel(MClust_FeatureNames{iC});
			ylabel(MClust_FeatureNames{jC});
		end
	end
end

return

function DrawAxes(xdim, ydim)

% -- get variables
global MClust_Clusters MClust_UnaccountedForOnly MClust_Colors MClust_Hide MClust_ClusterIndex MClust_FeatureData MClust_FDfn MClust_CurrentFeatures
global MClust_FeatureTimestamps MClust_FeatureNames

nClust = length(MClust_Clusters);
marker = '.';
markerSize = 1;
 
CurrentPointer = getptr(gcf); setptr(gcf, 'watch');
if ~strcmp(MClust_FeatureNames(xdim), MClust_CurrentFeatures(1))
   if strcmpi(MClust_FeatureNames{xdim}(1:4), 'time')
       MClust_FeatureData(:,1) = MClust_FeatureTimestamps;
       MClust_CurrentFeatures(1) = MClust_FeatureNames(xdim);
   else
       [fpath fname fext] = fileparts(MClust_FDfn);
       FeatureToGet = MClust_FeatureNames{xdim};
       FindColon = find(FeatureToGet == ':');
       temp = load(fullfile(fpath, [fname '_' FeatureToGet(1:FindColon-1) '.fd']),'-mat');
       FeatureIndex = strmatch(FeatureToGet,temp.FeatureNames);
       MClust_FeatureData(:,1) = temp.FeatureData(:,FeatureIndex);
       MClust_CurrentFeatures(1) = temp.FeatureNames(FeatureIndex); %MClust_ylbls(ydim);
   end;
end;
if ~strcmp(MClust_FeatureNames(ydim), MClust_CurrentFeatures(2))
   if strcmpi(MClust_FeatureNames{ydim}(1:4), 'time')
       MClust_FeatureData(:,2) = MClust_FeatureTimestamps;
       MClust_CurrentFeatures(2) = MClust_FeatureNames(ydim);
   else
       [fpath fname fext] = fileparts(MClust_FDfn);
       FeatureToGet = MClust_FeatureNames{ydim};
       FindColon = find(FeatureToGet == ':');
       temp = load(fullfile(fpath, [fname '_' FeatureToGet(1:FindColon-1) '.fd']),'-mat');
       FeatureIndex = strmatch(FeatureToGet,temp.FeatureNames);
       MClust_FeatureData(:,2) = temp.FeatureData(:,FeatureIndex);
       MClust_CurrentFeatures(2) = temp.FeatureNames(FeatureIndex); %MClust_ylbls(ydim);
   end;
end;
setptr(gcf, CurrentPointer{2});

for iC = 0:nClust     
    if ~MClust_Hide(iC+1)
        if iC == 0
            if MClust_UnaccountedForOnly
                MClust_ClusterIndex = ProcessClusters(MClust_FeatureData, MClust_Clusters);
                f = find(MClust_ClusterIndex == 0);                
                h = plot(MClust_FeatureData(f,1), MClust_FeatureData(f,2), marker);
            else
                h = plot(MClust_FeatureData(:,1), MClust_FeatureData(:,2), marker);
            end
        else         
            [f,MClust_Clusters{iC}] = FindInCluster(MClust_Clusters{iC}, MClust_FeatureData);
            h = plot(MClust_FeatureData(f,1), MClust_FeatureData(f,2), marker);
        end
        set(h, 'Color', MClust_Colors(iC+1,:));
        set(h, 'MarkerSize', markerSize);
    end
end
set(gca, 'XLim', [min(MClust_FeatureData(:,1)) max(MClust_FeatureData(:, 1))+0.0001], 'xtick', []);
set(gca, 'YLim', [min(MClust_FeatureData(:,2)) max(MClust_FeatureData(:, 2))+0.0001], 'ytick', []);
xlabel(MClust_CurrentFeatures{1});
ylabel(MClust_CurrentFeatures{2});

return