function [redraw, rekey, undoable] = FindBestCQProjections(iClust)

% [redraw, rekey, undoable] = FindBestCQProjections(iClust)
%
% Searches 2D projections to find the best axes for separating two clusters
% based on IsolationDistance and Lratio cluster quality measures.
% 
% If cluster 0 is selected as the second cluster, then the program attempts
% to separate the cluster from the noise.  Takes a long time if you have a
% lot of spikes.
%
% INPUTS
%    iClust
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder

% NCST 2003
% JCJ Oct 2004
% JCJ Feb 2008 Updated for MClust v3.5 to be a self-contained function with
%              callbacks handled by internal sub-function.
% JCJ April 2008 Reverted to using only features currently loaded in MClust
% JCJ April 2008 Further updated for MClust 3.5 compatibility



redraw   = false; % No changes are made to clusters
rekey    = false; % No changes are made to clusters
undoable = false; % No changes are made to clusters

if strcmpi(iClust,'Callback')
    cboHandle = gcbo;
    FBP_Callbacks(cboHandle);
    return
end


global MClust_FeatureTimestamps MClust_FeatureNames MClust_FeatureData MClust_FDfn MClust_Clusters

global MClust_FindBestSeparationFeatures MClust_Directory
% remove any blank features -- JCJ Sept 2007
MClust_FindBestSeparationFeatures=setdiff(MClust_FindBestSeparationFeatures,'');

if isempty(MClust_FindBestSeparationFeatures)
	for iF = 1:length(MClust_FeatureNames) - 1
		FeatureToGet = MClust_FeatureNames{iF};
		FindColon = find(FeatureToGet == ':');
		if isempty(strmatch(FeatureToGet(1:(FindColon - 1)),MClust_FindBestSeparationFeatures))
			MClust_FindBestSeparationFeatures{end+1} = FeatureToGet(1:(FindColon - 1));
		end
	end
end

featuresToUse = MClust_FindBestSeparationFeatures;

FindBestProjectionFigure = figure('Name','Find Best Projection','Tag','FindBestProjection','Units', 'Normalized');
%-------------------------------
% Alignment variables

uicHeight = 0.04;
uicWidth  = 0.25;
dX = 0.3;
XLocs = 0.1:dX:0.9;
dY = 0.04;
YLocs = 0.9:-dY:0.1;
FrameBorder = 0.01;

% Create Feature Listboxes
uicontrol('Parent', FindBestProjectionFigure,...
	'Style', 'text', 'String', 'FEATURES', 'Units', 'Normalized', 'Position', [XLocs(1) YLocs(4) 2*uicWidth uicHeight]);
uicontrol('Parent', FindBestProjectionFigure,...
	'Style', 'text', 'String', 'available', 'Units', 'Normalized', 'Position', [XLocs(1) YLocs(5) uicWidth uicHeight]);
ui_featuresIgnoreLB =  uicontrol('Parent', FindBestProjectionFigure,...
	'Units', 'Normalized', 'Position', [XLocs(1) YLocs(17) uicWidth 12*uicHeight],...
	'Style', 'listbox', 'Tag', 'FeaturesIgnoreListbox',...
	'Callback', 'TransferBetweenListboxes',... % Replaced 'FindBestCQProjections(''Callback'');MClustCallbacks',... April 2008 - JCJ - Simplifies functionality
	'HorizontalAlignment', 'right', ...
	'Enable','on', ...
	'TooltipString', 'These are features which are not included but are also available.');
uicontrol('Parent',FindBestProjectionFigure,...
	'Style', 'text', 'String', 'used', 'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth YLocs(5) uicWidth uicHeight]);
ui_featuresUseLB = uicontrol('Parent', FindBestProjectionFigure,...
	'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth YLocs(17) uicWidth 12*uicHeight],...
	'Style', 'listbox', 'Tag', 'FeaturesUseListbox',...
	'Callback', 'TransferBetweenListboxes',...% Replaced 'FindBestCQProjections(''Callback'');MClustCallbacks',... April 2008 - JCJ - Simplifies functionality
	'HorizontalAlignment', 'right', ...
	'Enable','on', ...
	'TooltipString', 'These features will be used for cluster separation.');
set(ui_featuresIgnoreLB, 'UserData', ui_featuresUseLB);
set(ui_featuresUseLB,    'UserData', ui_featuresIgnoreLB);

% Re-inserted April 2008 - JCJ - Simplifies functionality
featureFiles = get(findobj('Parent',findobj('Name','MClust'),'Tag', 'FeaturesUseListbox'),'String');

featureIgnoreString = {}; %featureFiles(end);
featureUseString = {};    %featureFiles(1:end-1);

for iF = 1:length(featureFiles)
	if strmatch(featureFiles{iF},MClust_FindBestSeparationFeatures)
		featureUseString{end+1} = featureFiles{iF};
	else
		featureIgnoreString{end+1} = featureFiles{iF};
	end
end

% Removed April 2008 - JCJ
% % % % Locate and load the names of feature files
% % % featureFiles =  sortcell(FindFiles('feature_*.m', ...
% % %     'StartingDirectory', fullfile(MClust_Directory, 'Features') ,'CheckSubdirs', 0));
% % % 
% % % featureIgnoreString= {};
% % % featureUseString = MClust_FindBestSeparationFeatures';
% % % 
% % % for iF = 1:length(featureFiles)
% % %    [dummy, featureFiles{iF}] = fileparts(featureFiles{iF});
% % %    featureFiles{iF} = featureFiles{iF}(9:end); % cut "feature_" off front for display
% % %    if isempty(strmatch(upper(featureFiles{iF}), upper(featureUseString),'exact'))
% % %       featureIgnoreString = cat(1, featureIgnoreString, featureFiles(iF));
% % %    end
% % % end

set(ui_featuresIgnoreLB, 'String', featureIgnoreString);
set(ui_featuresUseLB, 'String', featureUseString);

uicontrol('Parent',FindBestProjectionFigure, ...
	'Units', 'Normalized', 'Position', [XLocs(1) + uicWidth/2 YLocs(1) uicWidth/2 uicHeight], ...
	'Style', 'edit','Tag', 'Cluster1_FBS', 'String', num2str(iClust), ...
	'TooltipString', 'Cluster 1');	
uicontrol('Parent',FindBestProjectionFigure, ...
	'Units', 'Normalized', 'Position', [XLocs(1) + uicWidth/2 YLocs(2) uicWidth/2 uicHeight], ...
	'Style', 'Popupmenu','Tag', 'Cluster2_FBS', 'String', num2str((0:length(MClust_Clusters))'), ...
	'TooltipString', 'Cluster 2');

uicontrol('Parent',FindBestProjectionFigure, ...
	'Units', 'Normalized', 'Position', [XLocs(1) YLocs(20) uicWidth uicHeight], ...
    'Style', 'pushbutton','Tag', 'AcceptFeatures_FBS', 'String', 'Accept features', 'Callback', 'FindBestCQProjections(''Callback'')', ...
    'TooltipString', 'Accept features to use for finding the best 2D projection to separate these clusters');
uicontrol('Parent',FindBestProjectionFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(21) uicWidth uicHeight], ...
    'Style', 'pushbutton', 'String', 'Cancel', 'Callback', 'close');

% inserted April 2008 - JCJ - Simplifies functionality
set(findobj(gcf,'Tag', 'AcceptFeatures_FBS'),'enable','on')

uicontrol('Parent', FindBestProjectionFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(21) - uicHeight*11/10 uicWidth*3 uicHeight*8/10], ...
    'Style', 'text', 'String', 'Author: Jadin C. Jackson and Neil Schmitzer-Torbert;  jadincjackson@gmail.com', ...
    'TooltipString', 'Original version created by Neil Schmitzer-Torbert. Thanks Neil!');


% % % % Removed April 2008 - JCJ - Simplifies functionality
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % This forces the user to click on the Feature Selection box which triggers
% % % % MClustCallbacks to recalculate the MClust_Temp_FeatureNames varable... It
% % % % is a quick hack but it works. (Whithout it, the function uses only the
% % % % last features that the FindBestProjection, CutOnBestProjection,
% % % % RunBBclust, or RunKKwik functions used.)  -- JCJ Oct 2004
% % % set(findobj(gcf,'Tag', 'AcceptFeatures_FBS'),'enable','off')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FBP_Callbacks(cboHandle)

% FBP_Callbacks
%
% Searches 2D projections to find the best axes for separating two clusters
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
% JCJ Oct 2004


switch get(cboHandle, 'Tag')
    % % % % Removed April 2008 - JCJ - Simplifies functionality
    % % %     case {'FeaturesUseListbox';'FeaturesIgnoreListbox'}
    % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %         % This forces the user to click on the Feature Selection box which triggers
    % % %         % MClustCallbacks to recalculate the MClust_Temp_FeatureNames varable... It
    % % %         % is a quick hack but it works. (Whithout it, the function uses only the
    % % %         % last features that the FindBestProjection, CutOnBestProjection,
    % % %         % RunBBclust, or RunKKwik functions used.)  -- JCJ Oct 2004
    % % %         set(findobj(gcf,'Tag', 'AcceptFeatures_FBS'),'enable','on')
        
    case 'AcceptFeatures_FBS'
        global MClust_FindBestSeparationFeatures
        global MClust_FeatureTimestamps MClust_FeatureData MClust_FDfn MClust_Clusters
        
        nF=10;
        
        newFeatures = get(findobj('Parent',findobj('Name','Find Best Projection'),'Tag','FeaturesUseListbox','TooltipString', ...
            'These features will be used for cluster separation.'),'String');
        
%         iClust = str2num(get(findobj('TooltipString','Cluster 1'),'String'));
%         iClust2 = get(findobj('TooltipString','Cluster 2'),'Value');
        iClust = str2num(get(findobj('Tag', 'Cluster1_FBS'),'String'));
        iClust2 = get(findobj('Tag', 'Cluster2_FBS'),'Value');
        iClust2 = iClust2 - 1;
        if ~isempty(newFeatures)
            MClust_FindBestSeparationFeatures = {};
            MClust_FindBestSeparationFeatures = newFeatures;
            
            figure(findobj('Name','Find Best Projection'))
            close;
        else
            msgbox('Error: At least one feature must be selected');
        end
        
        FBS_Features = upper(MClust_FindBestSeparationFeatures);
        
        [f_1 MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});
        
        if iClust2 ~= 0
            [f_2 MClust_Clusters{iClust2}] = FindInCluster(MClust_Clusters{iClust2});
        else
            f_2 = find(~ismember(1:size(MClust_FeatureData,1),f_1))';
        end
        
        [Best_x, Best_y] =FBP(FBS_Features,f_1,f_2,nF);
        
        msg = {};
        msg{1} = ['Best Projections for: '];
        msg{2} = ['Cluster ' num2str(iClust) ' vs. Cluster ' num2str(iClust2)];
        msg{3} = ' ';
        for iP = 1:nF
            msg{end + 1} = [Best_x{iP} ' x ' Best_y{iP}];
        end
        msgbox(msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Best_x, Best_y] =FBP(FBS_Features,f_1,f_2,NBest)
% [Best_x, Best_y] =FBP(FBS_Features,f_1,f_2,NBest)
% 
%Searches all combinations of 2-D projections in FBS_Features to find the
%projection with the least overlap (assuming a Gaussian distribution) of
%cluster points with indices f_1 and f_2. Returns the NBest projections.
% 
%INPUTS
% FBS_Features    cell array of features names to find (e.g. {'energy'; 'peak'})
% f_1, f_2        indices into the two clusters to evaluate for overlap
% NBest           number of best projections to find
% 
%OUTPUTS
% Best_x, Best_y  Best x and y axis pairs
% 
% JCJ Sept 2004, adapted from code by NCST
% JCJ April 2008 updated for MClust 3.5 compatibility
 
global MClust_Temp_FeatureNames MClust_FeatureTimestamps MClust_FDfn 
global MClust_FDdn % JCJ added April 2008 for MClust 3.5 compatibility

if ~exist('MClust_Temp_FeatureNames','var') | isempty(MClust_Temp_FeatureNames);
    global MClust_FeatureNames
    MClust_Temp_FeatureNames=unique(MClust_FeatureNames);
end

MClust_Temp_FeatureNames=unique(MClust_Temp_FeatureNames);

nFeatures = length(MClust_Temp_FeatureNames) - 1;
mxDist = chi2inv(0.9999,2);
M_dists = repmat(nan,nFeatures);
FD_temp = repmat(nan,size(MClust_FeatureTimestamps,1),2);

FBS_Features=upper(FBS_Features);

if ~exist('NBest','var')
    NBest = nchoosek(length(FBS_Features),2);
end
    
[fpath fname fext] = fileparts(MClust_FDfn);
fpath = MClust_FDdn; % JCJ added April 2008 for MClust 3.5 compatibility

Curr_X = ' ';
Curr_Y = ' ';

pushdir(fpath);
% try
    DisplayProgress(0,nFeatures,'Title',['Finding Best Projections']);
    for iX = 1:nFeatures
        FeatureToGet = MClust_Temp_FeatureNames{iX};
        FindColon = find(FeatureToGet == ':');
        if strmatch(upper(FeatureToGet(1:FindColon-1)),FBS_Features,'exact')
            if ~strcmpi(FeatureToGet(1:FindColon-1),Curr_X)
                temp_X = load(fullfile(fpath, [fname '_' FeatureToGet(1:FindColon-1) '.fd']),'-mat');
                Curr_X = FeatureToGet(1:FindColon-1);
            end
            FeatureIndex = strmatch(FeatureToGet,temp_X.FeatureNames);
            FD_temp(:,1) = temp_X.FeatureData(:,FeatureIndex);
            for iY = (iX + 1):nFeatures
                FeatureToGet = MClust_Temp_FeatureNames{iY};
                FindColon = find(FeatureToGet == ':');
                if strmatch(upper(FeatureToGet(1:FindColon-1)),FBS_Features,'exact')
                    if ~strcmpi(FeatureToGet(1:FindColon-1),Curr_Y)
                        if ~exist([fname '_' FeatureToGet(1:FindColon-1) '.fd'],'file')
                            cd FD
                        end
                        temp_Y = load(fullfile(fpath, [fname '_' FeatureToGet(1:FindColon-1) '.fd']),'-mat');
                        Curr_Y = FeatureToGet(1:FindColon-1);
                    end
                    FeatureIndex = strmatch(FeatureToGet,temp_Y.FeatureNames);
                    FD_temp(:,2) = temp_Y.FeatureData(:,FeatureIndex);

                    Dists = mahal(FD_temp(f_2,:),FD_temp(f_1,:));

                    % Using a faster version of L
                    [H bins] = hist(log10(Dists),500);
                    Chi_CDF = 1 - chi2cdf(10.^bins,2);
                    M_dists(iX,iY) = -sum(H(:).*Chi_CDF(:));
                end
            end;
        end
        DisplayProgress(iX,nFeatures);
    end
    popdir;
% catch
%     popdir;
% end
Best_x = {}; Best_y = {};
for iP = 1:NBest
	[jnk yi] = max(max(M_dists));
	Best_y{iP} = MClust_Temp_FeatureNames{yi};
	[jnk xi] = max(M_dists(:,yi));
	Best_x{iP} = MClust_Temp_FeatureNames{xi};

	M_dists(xi,yi) = NaN;
end