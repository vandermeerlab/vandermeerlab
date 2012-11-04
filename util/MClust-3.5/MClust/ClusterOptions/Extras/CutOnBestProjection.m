function [redraw, rekey, undoable] = CutOnBestProjection(iClust)

% CutOnBestProjection(iClust)
%
% Finds a 2-D projection that is a superposition of the 3 best features to
% separate two clusters.
%
% If cluster 0 is selected as the second cluster, then the program attempts
% to separate the cluster from the noise.
%
% INPUTS
%    iClust
%
% OUTPUTS
%
% NONE
%
% CutOnBestProjection tries to find the 2-D projection that either 1)
% maximizes the distance between the clusters or 2) maximizes the distance
% while minimizing the overlap. The algorithm for 1) goes as follows:   
%
% * find the vector between the two cluster centers in the N-Dim feature
% space chosen by the user 
% * calculate the covariance matrix and it's Eigen vectors + values for
% each cluster 
% * Transform clusters into a coordinate system where cluster 1 is unit
% Gaussian centered at origin and covariance matrix of cluster 1 is
% diagonal (D1) (Basically what is done to get the Mahalanobis Distance).  
% -- this means there is a new vector from center of cluster1 to center of
% cluster2 
% -- covariance of cluster 2 in transformed into new coordinate system 
% * use the transformed covariance of cluster 2 to find the effective
% distance in each dimension from the cluster centers accounting for the
% covariances in cluster 2.  
% * then find the 3 dimensions that provide the largest separation between
% cluster centers given the distribution of cluster 2 
% * then since 3-D vector math is so easy find the cross product of the
% vector representing the largest axis of Cluster 2 with the vector between
% the two cluster centers to find a normal-vector which defines a plane
% slicing through the 3-D space.   
%
% Method 2) from above is accomplished similarly, however the 3-D space is
% immediately defined by the 3 top axes with minimum overlap (i.e. those
% with the smallest L-ratio) in 2-D projection combinations.  
% 
% If a user has FastICA on their Matlab path, then method 1) first performs
% an ICA before beginning. Method 2, finds the 6 top axes with minimum
% overlap (i.e. those with the smallest L-ratio) in 2-D projection
% combinations, then performs an ICA on these before finishing.  
%
%
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder

% NCST and JCJ Dec 2003
% JCJ Oct 2004
% JCJ Feb 2008 Updated for MClust v3.5 to be a self-contained function with
%              callbacks handled by internal sub-function.
% JCJ May 2010 Updated to ensure that the color of the pulldown font is not
%              the same as the color of the background

redraw    = false;
rekey     = false;
undoable  = false;

global MClust_FeatureTimestamps MClust_FeatureNames 
global MClust_FDfn MClust_Clusters
global CBP_data


if strcmpi(iClust,'Callback')
    cboHandle = gcbo;
    CBP_Callbacks(cboHandle);
    return
else
    if ~isa(MClust_Clusters{iClust},'mccluster')
        errordlg('Cluster must be of type ''mccluster'' to use CutOnBestProjection tool: First use ConvertCluster')
        return
    end
        
    % Search for already active CBP tool and abort if exist -- JCJ Oct 2007
    fh=findobj('Name','Cut On Best Projection');
    if ~isempty(fh)
        figure(fh);
        warndlg('Cut On Best Projection tool already in use.')
        return
    end
end


CBP_data = [];  % Important to clear every time since CBP_data is a global variable
CBP_data.iClust1 = [];
CBP_data.iClust2 = [];
CBP_data.frontback   = false;
CBP_data.Proj.Save   = false;
CBP_data.UseICA      = false; %Added ICA Sept 2007 JCJ
CBP_data.MaxUndoNmbr = 10;  %Added multi-step undo Oct 2007 JCJ
CBP_data.UseCQ       = false;

global MClust_CutOnBestProj MClust_Directory
if isempty(MClust_CutOnBestProj)
    for iF = 1:length(MClust_FeatureNames) - 1
        FeatureToGet = MClust_FeatureNames{iF};
        FindColon = find(FeatureToGet == ':');
        if isempty(strmatch(upper(FeatureToGet(1:(FindColon - 1))),upper(MClust_CutOnBestProj),'exact'))
            MClust_CutOnBestProj{end+1} = FeatureToGet(1:(FindColon - 1));
        end
    end
end

% remove any blank features
MClust_CutOnBestProj=setdiff(MClust_CutOnBestProj,'');

% featuresToUse = MClust_CutOnBestProj;

ScSize = get(0, 'ScreenSize');
CutOnBestProjFigure = figure('Name','Cut On Best Projection','Tag','CutOnBestProjection',...
    'Position',[ ScSize(3)/40 ScSize(4)*0.2 ScSize(3)*0.95 ScSize(4)*0.6]);

CBP_data.cbpfh=CutOnBestProjFigure;

set(CutOnBestProjFigure,'Units', 'Normalized');
%-------------------------------
% Alignment variables

uicHeight = 0.04;
uicWidth  = 0.25;
dX = uicWidth+uicWidth/5;
XLocs = 0.05:dX:0.95;
dY = 0.04;
YLocs = 0.9:-dY:0.1;
FrameBorder = 0.01;

%Added Stats Displays Sept 2007 JCJ
CuttingPlotPosition=[XLocs(2) YLocs(end) (XLocs(end)-XLocs(2))*2/3 YLocs(1)-YLocs(end)+dY];
XCorrPosition=[XLocs(2)+(XLocs(end)-XLocs(2))*13/18    YLocs(end)                           (XLocs(end)-XLocs(2))*1/3 (YLocs(1)-YLocs(end))*4/15];
HistISI2Position=[XLocs(2)+(XLocs(end)-XLocs(2))*13/18 YLocs(end)+(YLocs(1)-YLocs(end))*1/3+dY (XLocs(end)-XLocs(2))*1/3 (YLocs(1)-YLocs(end))*4/15];
HistISI1Position=[XLocs(2)+(XLocs(end)-XLocs(2))*13/18 YLocs(end)+(YLocs(1)-YLocs(end))*2/3+2*dY (XLocs(end)-XLocs(2))*1/3 (YLocs(1)-YLocs(end))*4/15];

%Added Stats Displays Sept 2007 JCJ
CBP_data.cbpfp.CuttingPlotPosition = CuttingPlotPosition;
CBP_data.cbpfp.XCorrPosition  = XCorrPosition;
CBP_data.cbpfp.HistISI1Position = HistISI1Position;
CBP_data.cbpfp.HistISI2Position = HistISI2Position;

% Create Feature Listboxes
uicontrol('Parent', CutOnBestProjFigure,...
    'Style', 'text', 'String', 'FEATURES', 'Units', 'Normalized', 'Position', [XLocs(1) YLocs(4) uicWidth uicHeight]);
uicontrol('Parent', CutOnBestProjFigure,...
    'Style', 'text', 'String', 'available', 'Units', 'Normalized', 'Position', [XLocs(1) YLocs(5) uicWidth/2 uicHeight]);
ui_featuresIgnoreLB =  uicontrol('Parent', CutOnBestProjFigure,...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(17) uicWidth/2 12*uicHeight],...
    'Style', 'listbox', 'Tag', 'FeaturesIgnoreListbox',...
    'Callback', 'CutOnBestProjection(''Callback'')',... %JCJ April 2008 Removed MClustCallbacks
    'HorizontalAlignment', 'right', ...
    'Enable','on', ...
    'TooltipString', 'These are features which are not included but are also available.');
uicontrol('Parent',CutOnBestProjFigure,...
    'Style', 'text', 'String', 'used', 'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth/2 YLocs(5) uicWidth/2 uicHeight]);
ui_featuresUseLB = uicontrol('Parent', CutOnBestProjFigure,...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth/2 YLocs(17) uicWidth/2 12*uicHeight],...
    'Style', 'listbox', 'Tag', 'FeaturesUseListbox',...
    'Callback', ';CutOnBestProjection(''Callback'')',...  %JCJ April 2008 Removed MClustCallbacks
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
    if strmatch(featureFiles{iF},MClust_CutOnBestProj)
        featureUseString{end+1} = featureFiles{iF};
    else
        featureIgnoreString{end+1} = featureFiles{iF};
    end
end

% Removed April 2008 - JCJ
% % %
% % % % Locate and load the names of feature files
% % % featureFiles =  sortcell(FindFiles('feature_*.m', ...
% % %     'StartingDirectory', fullfile(MClust_Directory, 'Features') ,'CheckSubdirs', 0));
% % %
% % % featureIgnoreString= {};
% % % featureUseString = MClust_CutOnBestProj';
% % %
% % % for iF = 1:length(featureFiles)
% % %     [dummy, featureFiles{iF}] = fileparts(featureFiles{iF});
% % %     featureFiles{iF} = featureFiles{iF}(9:end); % cut "feature_" off front for display
% % %     if isempty(strmatch(upper(featureFiles{iF}), upper(featureUseString),'exact'))
% % %         featureIgnoreString = cat(1, featureIgnoreString, featureFiles(iF));
% % %     end
% % % end

set(ui_featuresIgnoreLB, 'String', featureIgnoreString);
set(ui_featuresUseLB, 'String', featureUseString);

% Create cluster selection pulldowns
CBP_data.C1h = uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth/5  YLocs(1) uicWidth/6 uicHeight], ...
    'Style', 'Popupmenu','Tag', 'Cluster1pulldownCBP', 'String', num2str((1:length(MClust_Clusters))'), ...
    'TooltipString', 'Cluster A','Callback', 'CutOnBestProjection(''Callback'')');

set(CBP_data.C1h,'Value',iClust);
CBP_data.C2h=uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth/5 YLocs(2) uicWidth/6 uicHeight], ...
    'Style', 'Popupmenu','Tag', 'Cluster2pulldownCBP', 'String', num2str((0:length(MClust_Clusters))'), ...
    'TooltipString', 'Cluster B','Callback', 'CutOnBestProjection(''Callback'')');

% Create cutter buttons
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth*2/5 YLocs(1) + uicHeight*1/10 uicWidth/4 uicHeight*9/10], ...
    'Style', 'pushbutton','Tag', 'AddLimit', 'String', 'Add Limit', 'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'Add limit to Cluster A');
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth*2/5 YLocs(2) + uicHeight*1/10 uicWidth/4 uicHeight*9/10], ...
    'Style', 'pushbutton','Tag', 'AddLimit', 'String', 'Add Limit', 'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'Add limit to Cluster B');

set(findobj(gcf,'TooltipString', 'Add limit to Cluster A'),'enable','off')
set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')

uicontrol('Parent', CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(1) + uicHeight*1/10 uicWidth/5 uicHeight*9/10], ...
    'Style', 'checkbox','Value', 0, 'Tag', 'Hide1', 'String', 'Hide', ...
    'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'If checked, hides cluster specified in box 1.');
uicontrol('Parent', CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(2) + uicHeight*1/10 uicWidth/5 uicHeight*9/10], ...
    'Style', 'checkbox','Value', 0, 'Tag', 'Hide2', 'String', 'Hide', ...
    'Callback','CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'If checked, hides cluster specified in box 2.');

uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(2)-uicWidth/2 YLocs(1) + uicHeight*1/10 uicWidth*2/7 uicHeight*9/10], ...
    'Style', 'pushbutton','Tag', 'FrontBack', 'String', 'Front/Back', 'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'Toggles which cluster is in front.');
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(2)-uicWidth/2 YLocs(2) + uicHeight*1/10 uicWidth*2/7 uicHeight*9/10], ...
    'Style', 'pushbutton','Tag', 'Swap', 'String', 'Swap', 'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'Switches which cluster is being optimized with respect to the other.');
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(2)-uicWidth/2 YLocs(3) + uicHeight*1/10 uicWidth*2/7 uicHeight*9/10], ...
    'Style', 'checkbox','Value', 0, 'Tag', 'SaveProj', 'String', 'Save Proj', ...
    'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'If checked, keeps projection information.');

% Enable the ICA button if fastica function is found on path
ICAloc=which('fastica'); %Added Sept 2007 JCJ (+ Changes to Projection2DSuperposition and CutOnBestProjection(''Callback''))
if ~isempty(ICAloc)
    uicontrol('Parent',CutOnBestProjFigure, ...
        'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth*2/5 YLocs(3) + uicHeight*1/10 uicWidth*2/7 uicHeight*9/10], ...
        'Style', 'checkbox','Value', 0, 'Tag', 'ICA', 'String', 'Use ICA', ...
        'Callback', 'CutOnBestProjection(''Callback'')', ...
        'TooltipString', 'If checked, uses FastICA (independent component analysis) to pre-process feature space.');
else
    uicontrol('Parent',CutOnBestProjFigure, ...
        'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth*2/5 YLocs(3) + uicHeight*1/10 uicWidth*2/7 uicHeight*9/10], ...
        'Style', 'checkbox','Value', 0, 'Tag', 'ICA', 'String', 'Use ICA', ...
        'Callback', 'warndlg(''FastICA is required for the ICA function. Downloadable at http://www.cis.hut.fi/projects/ica'')', 'enable','off',...
        'TooltipString', 'FastICA is not found on the path. Download at http://www.cis.hut.fi/projects/ica');
end

uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(19) uicWidth/2 uicHeight], ...
    'Style', 'pushbutton','Tag', 'CalcBestProj', 'String', 'Recalc Proj', 'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'Find the best 2D projection to separate these clusters');
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(20) uicWidth/2 uicHeight], ...
    'Style', 'pushbutton', 'String', 'Accept Clusters','Tag','AcceptClusters', 'Callback', 'CutOnBestProjection(''Callback'')');
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(21) uicWidth/2 uicHeight], ...
    'Style', 'pushbutton', 'String', 'Cancel', 'Callback', 'close;clear global CBP_data');
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth/2 YLocs(19) uicWidth/4 uicHeight], ...
    'Style', 'pushbutton', 'String', 'Undo','Tag','Undo', 'Callback', 'CutOnBestProjection(''Callback'')',...
    'TooltipString', 'Steps backward, removing limits');
uicontrol('Parent',CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth*3/4 YLocs(19) uicWidth/4 uicHeight], ...
    'Style', 'pushbutton', 'String', 'Redo','Tag','Redo', 'Callback', 'CutOnBestProjection(''Callback'')',...
    'TooltipString', 'Steps forward, re-applying limits');

uicontrol('Parent', CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth*17/32 YLocs(20) + uicHeight*1/10 uicWidth*14/32 uicHeight*8/10], ...
    'Style', 'checkbox','Value', 1, 'Tag', 'Mdist', 'String', 'Separation', ...
    'Callback', 'CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'Uses Mahalanobis Distance to select dimensions');
uicontrol('Parent', CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1)+uicWidth*17/32 YLocs(21) + uicHeight*1/10 uicWidth*14/32 uicHeight*8/10], ...
    'Style', 'checkbox','Value', 0, 'Tag', 'CQ', 'String', 'Overlap', ...
    'Callback','CutOnBestProjection(''Callback'')', ...
    'TooltipString', 'Uses overlap of individual 2-D projections to select dimensions');

% JCJ April 2008 Update cluster data
CBP_Callbacks(CBP_data.C1h,'initialize')

uicontrol('Parent', CutOnBestProjFigure, ...
    'Units', 'Normalized', 'Position', [XLocs(1) YLocs(21) - uicHeight*11/10 uicWidth uicHeight*8/10], ...
    'Style', 'text', 'String', 'Author: Jadin C. Jackson;  jadincjackson@gmail.com', ...
    'TooltipString', 'includes code borrowed from Neil Schmitzer-Torbert. Thanks Neil!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CBP_Callbacks(cboHandle,SwitchTag)

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
% V3.0 JCJ Oct 2004
% V4.0 JCJ 2007 -- included cluster stats displays, undo, and ICA option
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data
global MClust_CutOnBestProj MClust_Colors
global MClust_FeatureTimestamps  MClust_FDfn MClust_Clusters
global MClust_Clusters MClust_TTdn MClust_TTfn MClust_TText MClust_ChannelValidity


StoreMClustUndoState = false;

Show1Handle = findobj(CBP_data.cbpfh, 'Tag', 'Hide1');
Show1Flag = ~get(Show1Handle, 'Value');
Show2Handle = findobj(CBP_data.cbpfh, 'Tag', 'Hide2');
Show2Flag = ~get(Show2Handle, 'Value');

%Added Stats Displays Sept 2007 JCJ
CuttingPlotPosition=CBP_data.cbpfp.CuttingPlotPosition;
XCorrPosition=CBP_data.cbpfp.XCorrPosition;
HistISI1Position=CBP_data.cbpfp.HistISI1Position;
HistISI2Position=CBP_data.cbpfp.HistISI2Position;

if nargin<2
    SwitchTag=get(cboHandle, 'Tag');
end

%%% Close any open warning dialogs
openFigures=findobj(allchild(0),'flat','Visible','on');
Names=get(openFigures,'Name');
close(openFigures(strncmpi(Names,'CBP',3)))

switch SwitchTag
    case{'FeaturesIgnoreListbox';'FeaturesUseListbox'}
        CBP_data.Proj.Save = false;
        TransferBetweenListboxes;  %Added April 2008 JCJ - simplifies functionality, adds speed in MClust 3.5
        set(findobj(gcf,'Tag', 'SaveProj'),'Value',0);
    case 'SaveProj'
        CBP_data.Proj.Save = get(findobj(gcf,'Tag', 'SaveProj'),'Value');
    case 'Mdist'
        if get(findobj(gcf,'Tag', 'Mdist'),'Value')
            set(findobj(gcf,'Tag', 'CQ'),'Value',0)
            CBP_data.UseCQ = false;
        else
            set(findobj(gcf,'Tag', 'CQ'),'Value',1)
            CBP_data.UseCQ = true;
        end
    case 'CQ'
        if get(findobj(gcf,'Tag', 'CQ'),'Value')
            set(findobj(gcf,'Tag', 'Mdist'),'Value',0)
            CBP_data.UseCQ = true;
        else
            set(findobj(gcf,'Tag', 'Mdist'),'Value',1)
            CBP_data.UseCQ = false;
        end
    case 'ICA'  %Added Sept 2007 JCJ (+ Changes to CutOnBestProjection and Projection2DSuperposition)

        CBP_data.UseICA = get(findobj(gcf,'Tag', 'ICA'),'Value');

    case 'Hide1'

        if Show1Flag
            set(CBP_data.h(1),'Visible','on')
        else
            set(CBP_data.h(1),'Visible','off')
        end
        %     axis(CBP_data.a)

    case 'Hide2'
        if Show2Flag
            set(CBP_data.h(2),'Visible','on')
        else
            set(CBP_data.h(2),'Visible','off')
        end
        %     axis(CBP_data.a)

    case 'FrontBack'

        global MClust_Colors

        [Marker, MarkerSize] = GetMarker;

        figure(findobj('Name','Cut On Best Projection'))
        subplot('Position',CuttingPlotPosition)
        if ~CBP_data.frontback
            delete(CBP_data.h(2))
            CBP_data.h(2)=plot(CBP_data.x2(CBP_data.keeps2),CBP_data.y2(CBP_data.keeps2),...
                Marker,'MarkerSize',MarkerSize,...
                'Color',MClust_Colors(CBP_data.iClust2 + 1,:));
        else
            delete(CBP_data.h(1))
            CBP_data.h(1)=plot(CBP_data.x1(CBP_data.keeps1),CBP_data.y1(CBP_data.keeps1),...
                Marker,'MarkerSize',MarkerSize,...
                'Color',MClust_Colors(CBP_data.iClust1 + 1,:));
        end
        axis(CBP_data.a)
        CBP_data.frontback=~CBP_data.frontback;
    case 'CalcBestProj'
        % Turn off buttons temporarily to reduce conflicts during calculation -- JCJ Sept 2007
        Sub_DisableButtons
        set(findobj(gcf,'TooltipString', 'Add limit to Cluster A'),'enable','off')
        set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
        set(findobj(gcf,'Tag', 'Swap'),'enable','off')

        if all(isfield(CBP_data,{'a','h','x1','y1','x2','y2','keeps1','keeps2','Proj'}))
            Sub_StoreUNDO;
        end

        newFeatures = get(findobj('Parent',findobj('Name','Cut On Best Projection'),'Tag','FeaturesUseListbox','TooltipString', ...
            'These features will be used for cluster separation.'),'String');

        %     iClust = get(findobj('TooltipString','Cluster A'),'Value'); %str2num(get(findobj('TooltipString','Cluster A'),'String'));
        %     iClust2 = get(findobj('TooltipString','Cluster B'),'Value');
        iClust = get(findobj('Tag', 'Cluster1pulldownCBP'),'Value'); %str2num(get(findobj('TooltipString','Cluster A'),'String'));
        iClust2 = get(findobj('Tag', 'Cluster2pulldownCBP'),'Value');
        %     iClust = iClust - 1;
        iClust2 = iClust2 - 1;

        [f1 MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});

        if iClust2 ~= 0
            [f2 MClust_Clusters{iClust2}] = FindInCluster(MClust_Clusters{iClust2});
        else
            f2 = find(~ismember(1:size(MClust_FeatureTimestamps,1),f1))';
            set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
            set(findobj(gcf,'Tag', 'Swap'),'enable','off')

        end

        if ~isempty(newFeatures)
            MClust_CutOnBestProj = {};
            MClust_CutOnBestProj = newFeatures;
            % 	close;
        else
            msgbox('Error: At least one feature must be selected');
        end

        figure(findobj('Name','Cut On Best Projection'))
        %     subplot(1,20,7:17);
        subplot('Position',CuttingPlotPosition)
        Projection2dSuperposition(f1,f2,[MClust_TTdn filesep MClust_TTfn MClust_TText],MClust_ChannelValidity,iClust,iClust2,'FeaturesToGet',MClust_CutOnBestProj,'CuttingPlotPosition',CuttingPlotPosition);

        if Show1Flag
            set(CBP_data.h(1),'Visible','on')
        else
            set(CBP_data.h(1),'Visible','off')
        end

        if Show2Flag
            set(CBP_data.h(2),'Visible','on')
        else
            set(CBP_data.h(2),'Visible','off')
        end
        axis(CBP_data.a)
        Sub_UpdateStats(XCorrPosition,HistISI1Position,HistISI2Position) %Added Stats Displays Sept 2007 JCJ

        % Moved to reduce conflicts during calculation -- JCJ Sept 2007
        set(findobj(gcf,'TooltipString', 'Add limit to Cluster A'),'enable','on')
        if iClust2 ~= 0
            set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','on')
            set(findobj(gcf,'Tag', 'Swap'),'enable','on')
        end
        %Re-enable other operations -- Sept 2007 JCJ
        Sub_EnableButtons

    case 'AddLimit'


        global MClust_Colors

        [Marker, MarkerSize] = GetMarker;

        Sub_StoreUNDO

        %Add Limit
        LimitError=true;
        while LimitError
            try
                figure(findobj('Name','Cut On Best Projection'))
                subplot('Position',CuttingPlotPosition);
                [x y] = DrawConvexHull;
                LimitError=false;
            end
        end

        % Limit limits to within boundary of axes ~JCJ Oct 2007
        a=axis;
        x=max(x,a(1));
        x=min(x,a(2));
        y=max(y,a(3));
        y=min(y,a(4));

        if strmatch('Add limit to Cluster A',get(cboHandle,'Tooltipstring'))
            iClust = CBP_data.iClust1;
            dx = CBP_data.x1;
            dy = CBP_data.y1;
            keeps = CBP_data.keeps1;
        else
            iClust = CBP_data.iClust2;
            dx = CBP_data.x2;
            dy = CBP_data.y2;
            keeps = CBP_data.keeps2;
        end

        if length(x > 2)
            f0 = find(~InPolygon(dx, dy, x, y) & keeps);
            %         MClust_Clusters{iClust} = Remove_Points(MClust_Clusters{iClust},f1(f0));

            plot(x,y,'Color',MClust_Colors(iClust + 1,:),'LineWidth',2);
            plot(dx(f0),dy(f0),'*k','MarkerSize',MarkerSize);
        end

        CBP_data.a=axis;

        if strmatch('Add limit to Cluster A',get(cboHandle,'Tooltipstring'))
            CBP_data.keeps1(f0) = 0;
        else
            CBP_data.keeps2(f0) = 0;
        end

        Sub_UpdateStats(XCorrPosition,HistISI1Position,HistISI2Position) %Added Stats Displays Sept 2007 JCJ



    case {'Undo';'Redo'}
        global MClust_Colors

        [Marker, MarkerSize] = GetMarker;
        Sub_DisableButtons
        Sub_ClearPlots(CuttingPlotPosition,XCorrPosition,HistISI1Position,HistISI2Position)

        if strcmp(SwitchTag,'Undo')
            Sub_RecallUndo;
        elseif strcmp(SwitchTag,'Redo')
            Sub_REDO;
        else
            error('Corrupted input')
        end

        f1 = CBP_data.f1(find(CBP_data.keeps1));
        f2 = CBP_data.f2(find(CBP_data.keeps2));

        idxS=ismember(f1,f2);

        figure(findobj('Name','Cut On Best Projection'))
        subplot('Position',CuttingPlotPosition)
        if CBP_data.frontback
            CBP_data.h(1)=plot(CBP_data.x1(CBP_data.keeps1),CBP_data.y1(CBP_data.keeps1),...
                Marker,'MarkerSize',MarkerSize,...
                'Color',MClust_Colors(CBP_data.iClust1 + 1,:));
            hold on
            CBP_data.h(2)=plot(CBP_data.x2(CBP_data.keeps2),CBP_data.y2(CBP_data.keeps2),...
                Marker,'MarkerSize',MarkerSize,...
                'Color',MClust_Colors(CBP_data.iClust2 + 1,:));
        else
            CBP_data.h(2)=plot(CBP_data.x2(CBP_data.keeps2),CBP_data.y2(CBP_data.keeps2),...
                Marker,'MarkerSize',MarkerSize,...
                'Color',MClust_Colors(CBP_data.iClust2 + 1,:));
            hold on
            CBP_data.h(1)=plot(CBP_data.x1(CBP_data.keeps1),CBP_data.y1(CBP_data.keeps1),...
                Marker,'MarkerSize',MarkerSize,...
                'Color',MClust_Colors(CBP_data.iClust1 + 1,:));
        end
        plot(CBP_data.x1(idxS),CBP_data.y1(idxS),'og','MarkerSize',MarkerSize+4);

        if Show1Flag
            set(CBP_data.h(1),'Visible','on')
        else
            set(CBP_data.h(1),'Visible','off')
        end

        if Show2Flag
            set(CBP_data.h(2),'Visible','on')
        else
            set(CBP_data.h(2),'Visible','off')
        end

        axis(CBP_data.a)
        Sub_EnableButtons

        Sub_UpdateStats(XCorrPosition,HistISI1Position,HistISI2Position) %Added Stats Displays Sept 2007 JCJ

    case 'AcceptClusters'
        global MClust_Clusters

        %         MClustCutterStoreUndo('Cut On Best Proj')
        MClustUndoState = MClustCutterUndoGetCurrentState(['Cluster ' num2str(CBP_data.iClust1) ' vs ' 'Cluster ' num2str(CBP_data.iClust2) ' : CutOnBestProj']);
        StoreMClustUndoState  = true;

        MClust_Clusters{CBP_data.iClust1} = Restrict_Points(MClust_Clusters{CBP_data.iClust1},CBP_data.f1(find(~CBP_data.keeps1)));
        if CBP_data.iClust2 ~= 0
            MClust_Clusters{CBP_data.iClust2} = Restrict_Points(MClust_Clusters{CBP_data.iClust2},CBP_data.f2(find(~CBP_data.keeps2)));
        end
        figure(findobj('Name','Cut On Best Projection'))
        close;

        CHDrawingAxisWindow = findobj('Type','figure','Tag', 'CHDrawingAxisWindow');
        if ~isempty(CHDrawingAxisWindow)
            close(CHDrawingAxisWindow);
        end
        MClustCutterCallbacks('RedrawAxes')


    case {'Swap'}  %Separated from pulldown update detection Sept 2007 JCJ
        CurriClust2 = get(findobj('Tag', 'Cluster2pulldownCBP'),'Value'); %str2num(get(findobj('TooltipString','Cluster A'),'String'));

        if CurriClust2 ~= 0
            %Temporarily disable other operations -- Sept 2007 JCJ
            Sub_DisableButtons

            % Swap cluster numbers
            iClust2 = get(findobj('Tag', 'Cluster1pulldownCBP'),'Value'); %str2num(get(findobj('TooltipString','Cluster A'),'String'));
            iClust1 = get(findobj('Tag', 'Cluster2pulldownCBP'),'Value');
            %     iClust1 = iClust1 - 1;
            iClust1 = iClust1 - 1;

            CBP_data.iClust1=iClust1;
            CBP_data.iClust2=iClust2;

            % Swap cluster attributes
            Temp=CBP_data;

            CBP_data.x1=Temp.x2;
            CBP_data.y1=Temp.y2;
            CBP_data.f1=Temp.f2;
            CBP_data.keeps1=Temp.keeps2;

            CBP_data.x2=Temp.x1;
            CBP_data.y2=Temp.y1;
            CBP_data.f2=Temp.f1;
            CBP_data.keeps2=Temp.keeps1;

            % Clear UNDO (since hitting UNDO now would mess everything up)
            Sub_SwapUNDO

        else
            warndlg(['Cannot Swap when Cluster B is zero cluster.'], 'Swap unsuccessful.');
            set(findobj(gcf,'Tag', 'Swap'),'enable','off')
        end

        % Update Display and Menus
        Sub_ClearPlots(CuttingPlotPosition,XCorrPosition,HistISI1Position,HistISI2Position)
        Sub_UpdateStats(XCorrPosition,HistISI1Position,HistISI2Position) %Added Stats Displays Sept 2007 JCJ

        set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
        set(findobj(gcf,'TooltipString', 'Add limit to Cluster A'),'enable','off')

        set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',iClust1);
        set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',iClust2+1);

        % JCJ April 2008 - add color to pull down to help identify clusters with points and MClustCutter.
        % JCJ May 2010 - make sure that the color of the font is not the
        %                same as the color of the background
        Sub_UpdatePullDownColor;
        
        %Re-enable other operations -- Sept 2007 JCJ
        Sub_EnableButtons;

    case {'initialize'}
        % JCJ April 2008 - added to better initialize variables and display

        iClust1 = get(findobj('Tag', 'Cluster1pulldownCBP'),'Value'); %str2num(get(findobj('TooltipString','Cluster A'),'String'));
        iClust2 = get(findobj('Tag', 'Cluster2pulldownCBP'),'Value');
        iClust2 = iClust2 - 1;
        
        %%% Get data for cluster 1
        [f1 MClust_Clusters{iClust1}] = FindInCluster(MClust_Clusters{iClust1});

        CBP_data.iClust1 = iClust1;
        CBP_data.f1 = f1;
        CBP_data.keeps1 = logical(repmat(1,size(f1)));
        CBP_data.x1 = repmat(NaN,size(f1));
        CBP_data.y1 = repmat(NaN,size(f1));


        %%% Get data for cluster 2
        if iClust2 ~= 0 % added to stop error if user sets cluster B back to zero - JCJ Sept 2007
            set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','on')
            set(findobj(gcf,'Tag', 'Swap'),'enable','on')

            [f2 MClust_Clusters{iClust2}] = FindInCluster(MClust_Clusters{iClust2});
        else
            f2 = find(~ismember(1:size(MClust_FeatureTimestamps,1),CBP_data.f1))';
            set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
            set(findobj(gcf,'Tag', 'Swap'),'enable','off')

        end

        CBP_data.iClust2 = iClust2;
        CBP_data.f2 = f2;
        CBP_data.keeps2 = logical(repmat(1,size(f2)));
        CBP_data.x2 = repmat(NaN,size(f2));
        CBP_data.y2 = repmat(NaN,size(f2));

        %     subplot(1,20,7:17);
        Sub_EnableButtons;

        Sub_ClearPlots(CuttingPlotPosition,XCorrPosition,HistISI1Position,HistISI2Position)
        set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
        set(findobj(gcf,'TooltipString', 'Add limit to Cluster A'),'enable','off')
        set(Show1Handle,'enable','off')
        set(Show2Handle,'enable','off')
        
        % JCJ May 2010 - make sure that the color of the font is not the
        %                same as the color of the background
        Sub_UpdatePullDownColor;

    case {'Cluster1pulldownCBP';'Cluster2pulldownCBP'}

        % Detect which cluster(s) changed
        iClust1 = get(findobj('Tag', 'Cluster1pulldownCBP'),'Value'); %str2num(get(findobj('TooltipString','Cluster A'),'String'));
        iClust2 = get(findobj('Tag', 'Cluster2pulldownCBP'),'Value');
        %     iClust1 = iClust1 - 1;
        iClust2 = iClust2 - 1;
        ClusterChange = false; % initialize -- if the clusters changed clear Undo - JCJ Oct 2007

        global MClust_Clusters 

        if isempty(CBP_data.iClust1)
            return
        end

        if iClust1 == iClust2
            %         set(findobj('TooltipString','Cluster A'),'Value',CBP_data.iClust1);
            %         set(findobj('TooltipString','Cluster B'),'Value',CBP_data.iClust2 + 1);
            set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',CBP_data.iClust1);
            set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',CBP_data.iClust2 + 1);
            return
        end

        newClust1 = ~(CBP_data.iClust1 == iClust1);
        newClust2 = ~(CBP_data.iClust2 == iClust2);

        if newClust1
            if ~isa(MClust_Clusters{iClust1},'mccluster')
                set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',CBP_data.iClust1);
                set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',CBP_data.iClust2 + 1);

                errordlg(['Cluster ' num2str(iClust1) ' must be of type ''mccluster'' to use CutOnBestProjection tool: First use ConvertCluster'])
                return
            end

            if ~isempty(strmatch('keeps1',fields(CBP_data)))
                if any(~CBP_data.keeps1)
                    ynClose = questdlg(['Accept changes to Cluster ' num2str(CBP_data.iClust1) '?.  No undo available. Are you sure?'], 'New Cluster A', 'Yes', 'No', 'Cancel', 'Cancel');
                    if strcmp(ynClose,'Cancel')
                        set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',CBP_data.iClust1);
                        set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',CBP_data.iClust2 + 1);
                        return
                    elseif strcmp(ynClose,'Yes')
                        MClustUndoState = MClustCutterUndoGetCurrentState(['Cluster ' num2str(CBP_data.iClust1) ' vs ' 'Cluster ' num2str(CBP_data.iClust2) ' : CutOnBestProj']);
                        StoreMClustUndoState  = true;

                        MClust_Clusters{CBP_data.iClust1} = Restrict_Points(MClust_Clusters{CBP_data.iClust1},CBP_data.f1(find(~CBP_data.keeps1)));
                        ClusterChange = true; % if the clusters changed clear Undo - JCJ Oct 2007

                    end
                end
                [f1 MClust_Clusters{iClust1}] = FindInCluster(MClust_Clusters{iClust1});

                CBP_data.iClust1 = iClust1;
                CBP_data.f1 = f1;
                CBP_data.keeps1 = logical(repmat(1,size(f1)));
                CBP_data.x1 = repmat(NaN,size(f1));
                CBP_data.y1 = repmat(NaN,size(f1));

                %     subplot(1,20,7:17);
                Sub_EnableButtons

                Sub_ClearPlots(CuttingPlotPosition,XCorrPosition,HistISI1Position,HistISI2Position)
                set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
                set(findobj(gcf,'TooltipString', 'Add limit to Cluster A'),'enable','off')
            end
        end

        if newClust2 
            if iClust2 ~= 0 && ~isa(MClust_Clusters{iClust2},'mccluster')
                set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',CBP_data.iClust1);
                set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',CBP_data.iClust2 + 1);

                errordlg(['Cluster ' num2str(iClust2) ' must be of type ''mccluster'' to use CutOnBestProjection tool: First use ConvertCluster'])
                return
            end

            if ~isempty(strmatch('keeps1',fields(CBP_data)))
                if any(~CBP_data.keeps2) & CBP_data.iClust2 ~= 0
                    ynClose = questdlg(['Accept changes to Cluster ' num2str(CBP_data.iClust2) '?.  No undo available. Are you sure?'], 'New Cluster B', 'Yes', 'No', 'Cancel', 'Cancel');
                    if strcmp(ynClose,'Cancel')
                        set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',CBP_data.iClust1);
                        set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',CBP_data.iClust2 + 1);
                        return
                    elseif strcmp(ynClose,'Yes')
                        MClustUndoState = MClustCutterUndoGetCurrentState(['Cluster ' num2str(CBP_data.iClust1) ' vs ' 'Cluster ' num2str(CBP_data.iClust2) ' : CutOnBestProj']);
                        StoreMClustUndoState  = true;

                        MClust_Clusters{CBP_data.iClust2} = Restrict_Points(MClust_Clusters{CBP_data.iClust2},CBP_data.f2(find(~CBP_data.keeps2)));
                        ClusterChange = true; % if the clusters changed clear Undo - JCJ Oct 2007
                    end
                end
                if iClust2 ~= 0 % added to stop error if user sets cluster B back to zero - JCJ Sept 2007
                    set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','on')
                    set(findobj(gcf,'Tag', 'Swap'),'enable','on')

                    [f2 MClust_Clusters{iClust2}] = FindInCluster(MClust_Clusters{iClust2});
                else
                    f2 = find(~ismember(1:size(MClust_FeatureTimestamps,1),CBP_data.f1))';
                    set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
                    set(findobj(gcf,'Tag', 'Swap'),'enable','off')

                end

                CBP_data.iClust2 = iClust2;
                CBP_data.f2 = f2;
                CBP_data.keeps2 = logical(repmat(1,size(f2)));
                CBP_data.x2 = NaN(size(f2));
                CBP_data.y2 = NaN(size(f2));

                %     subplot(1,20,7:17);
                Sub_ClearPlots(CuttingPlotPosition,XCorrPosition,HistISI1Position,HistISI2Position)
                Sub_EnableButtons
                set(findobj(gcf,'TooltipString', 'Add limit to Cluster B'),'enable','off')
                set(findobj(gcf,'TooltipString', 'Add limit to Cluster A'),'enable','off')
            end
        end
        
        %%%%%%%%%%%%%%%%%
        % JCJ April 2008 - add color to pull down to help identify clusters with points and MClustCutter.
        % JCJ May 2010 - make sure that the color of the font is not the
        %                same as the color of the background
        Sub_UpdatePullDownColor;

        % Update Display and Menus
        Sub_UpdateStats(XCorrPosition,HistISI1Position,HistISI2Position) %Added Stats Displays Sept 2007 JCJ

        if ClusterChange % if the clusters changed clear Undo - JCJ Oct 2007
            Sub_ClearUNDO;
        end
end

%%
if CBP_data.Proj.Save
    set(findobj(gcf,'Tag', 'Mdist'),'enable','off')
    set(findobj(gcf,'Tag', 'CQ'),'enable','off')
else
    set(findobj(gcf,'Tag', 'Mdist'),'enable','on')
    set(findobj(gcf,'Tag', 'CQ'),'enable','on')
end

%%% Check cluster list against MClust master list
if   ~strcmpi(SwitchTag,'AcceptClusters') && ((size(get(CBP_data.C1h,'String'),1)~=length(MClust_Clusters)) || (size(get(CBP_data.C2h,'String'),1)~=(length(MClust_Clusters)+1)))
    warndlg('Number of MClust clusters changed. This may cause update errors or data corruption if current cluster identities change in MClust.','CBP Warning: MClust Clusters Changed')
    set(CBP_data.C1h,'String', num2str((1:length(MClust_Clusters))'),'Value',CBP_data.iClust1);
    set(CBP_data.C2h,'String', num2str((0:length(MClust_Clusters))'),'Value', CBP_data.iClust2 + 1);
end

if  StoreMClustUndoState
    MClustCutterUndoStore(MClustUndoState);
end



function Sub_UpdatePullDownColor
% JCJ May 2010 - make sure that the color of the font is not the
%                same as the color of the background
%%
global CBP_data;
global MClust_Colors;

contrastthresh = 0.1;
BGcolor = MClust_Colors(CBP_data.iClust1+1,:);
FGcolor = 1-MClust_Colors(CBP_data.iClust1+1,:);
contrast=sqrt(BGcolor*FGcolor')/norm(FGcolor,2);

if abs(contrast-1)<contrastthresh
    % if the contrast is too low then set the font color to white
    % if the background is darker (less) than [0.5 0.5 0.5] or black if the
    % background is lighter (greater than).
    FGcolor = [1 1 1]*double(norm(BGcolor)<(sqrt(3)/2));
end
set(findobj('Tag', 'Cluster1pulldownCBP'),'BackgroundColor',BGcolor,'ForegroundColor',FGcolor);

BGcolor = MClust_Colors(CBP_data.iClust2+1,:);
FGcolor = 1-MClust_Colors(CBP_data.iClust2+1,:);
contrast=sqrt(BGcolor*FGcolor')/norm(FGcolor,2);

if abs(contrast-1)<contrastthresh
    % if the contrast is too low then set the font color to white
    % if the background is darker (less) than [0.5 0.5 0.5] or black if the
    % background is lighter (greater than).
    FGcolor = [1 1 1]*double(norm(BGcolor)<(sqrt(3)/2));
end
set(findobj('Tag', 'Cluster2pulldownCBP'),'BackgroundColor',BGcolor,'ForegroundColor',FGcolor);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Added Button Disable -- Sept 2007 JCJ
function Sub_DisableButtons
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data
%Disable other operations -- Sept 2007 JCJ
set(findobj(CBP_data.cbpfh,'Tag', 'CalcBestProj'),'enable','off')
set(findobj(CBP_data.cbpfh,'Tag', 'AcceptClusters'),'enable','off')
set(findobj(CBP_data.cbpfh,'Tag', 'SaveProj'),'enable','off')
set(findobj(CBP_data.cbpfh,'Tag', 'ICA'),'enable','off')
set(findobj(CBP_data.cbpfh,'Tag', 'Undo'),'enable','off')
set(findobj(CBP_data.cbpfh,'Tag', 'Redo'),'enable','off')
set(findobj(CBP_data.cbpfh,'Tag', 'CQ'),'enable','off')
set(findobj(CBP_data.cbpfh,'Tag', 'Mdist'),'enable','off')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Added Button Enable -- Sept 2007 JCJ
function Sub_EnableButtons
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data
%Re-enable other operations -- Sept 2007 JCJ
set(findobj(CBP_data.cbpfh,'Tag', 'CalcBestProj'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'AcceptClusters'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'SaveProj'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'ICA'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'Undo'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'Redo'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'CQ'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'Mdist'),'enable','on')
set(findobj(CBP_data.cbpfh,'Tag', 'Hide1'),'enable','on')  % JCJ April 2008
set(findobj(CBP_data.cbpfh,'Tag', 'Hide2'),'enable','on')  % JCJ April 2008



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Added Stats Displays Sept 2007 JCJ
function Sub_ClearPlots(CuttingPlotPosition,XCorrPosition,HistISI1Position,HistISI2Position)
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data
figure(CBP_data.cbpfh)
subplot('Position',CuttingPlotPosition);hold off;plot(nan,nan)
subplot('Position',XCorrPosition);hold off;plot(nan,nan)
subplot('Position',HistISI1Position);hold off;plot(nan,nan)
subplot('Position',HistISI2Position);hold off;plot(nan,nan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Added Stats Displays Sept 2007 JCJ
function Sub_UpdateStats(XCorrPosition,HistISI1Position,HistISI2Position)
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
threshold=2; %2 ms
binsize_msec=1;
xcorr_width_msec=500;

global CBP_data MClust_TTData MClust_Clusters  MClust_FeatureTimestamps

f1=FindInCluster(MClust_Clusters{CBP_data.iClust1});
if CBP_data.iClust2
    f2=FindInCluster(MClust_Clusters{CBP_data.iClust2});
end

f1n=FindInCluster(Restrict_Points(MClust_Clusters{CBP_data.iClust1},CBP_data.f1(find(~CBP_data.keeps1))));
if CBP_data.iClust2
    f2n=FindInCluster(Restrict_Points(MClust_Clusters{CBP_data.iClust2},CBP_data.f2(find(~CBP_data.keeps2))));
end

%%%%%%%%%%%%%X-corr%%%%%%%%%%%%%%%%%
figure(CBP_data.cbpfh)
subplot('Position',XCorrPosition);hold off;plot(nan,nan)
if CBP_data.iClust2
    ts1 = Range(ExtractCluster(MClust_TTData, f1),'ts');
    ts2 = Range(ExtractCluster(MClust_TTData, f2),'ts');
    ts1n = Range(ExtractCluster(MClust_TTData, f1n),'ts');
    ts2n = Range(ExtractCluster(MClust_TTData, f2n),'ts');
    nbins = round(xcorr_width_msec/binsize_msec);
    figure(CBP_data.cbpfh)
    subplot('Position',XCorrPosition);
    if isempty(ts1) | isempty(ts2)
        msgbox('No spikes in one of the chosen clusters.')
    else
        [Y, xdim] = CrossCorr(ts1, ts2, binsize_msec, nbins);
        plot(xdim,Y,'r')
        hold on
        [Y, xdim] = CrossCorr(ts1n, ts2n, binsize_msec, nbins);
        plot(xdim,Y)
        axis tight
        xlabel('msec')
        title(['XCorr ' num2str(CBP_data.iClust1) ' Vs ' num2str(CBP_data.iClust2) ]);
        legend('old','new');legend('boxoff')
    end
    zoom on
else
    title('No data displayed when using cluster 0')
end

%%%%%%%%%%%%%Hist ISI Cluster A%%%%%%%%%%%%%%%%%
figure(CBP_data.cbpfh)
subplot('Position',HistISI1Position);hold off;plot(nan,nan)
HistISI(ts(MClust_FeatureTimestamps(f1)));c=get(gca,'Children');set(c(2),'Color','r');hold on
HistISI(ts(MClust_FeatureTimestamps(f1n)));
axis tight
a=axis;
hold on
hold on;plot([1 1]*threshold,[a(3:4)],'g:')
title(['Cluster ' num2str(CBP_data.iClust1) ' ISI Histogram'])% num2str(sum(H(N<threshold))+1) '/' num2str(sum(H)) ' ISIs less than ' num2str(threshold) 'ms']);
xlabel([]);
zoom on

%%%%%%%%%%%%%Hist ISI Cluster B%%%%%%%%%%%%%%%%%
figure(CBP_data.cbpfh)
subplot('Position',HistISI2Position);hold off;plot(nan,nan)
if CBP_data.iClust2
    HistISI(ts(MClust_FeatureTimestamps(f2)));c=get(gca,'Children');set(c(2),'Color','r');hold on
    HistISI(ts(MClust_FeatureTimestamps(f2n)));
    axis tight
    a=axis;
    hold on
    hold on;plot([1 1]*threshold,[a(3:4)],'g:')
    title(['Cluster ' num2str(CBP_data.iClust2) ' ISI Histogram'])%  ' num2str(sum(H(N<threshold))+1) '/' num2str(sum(H)) ' ISIs less than ' num2str(threshold) 'ms']);
    xlabel([]);
    zoom on
else
    title('No data displayed when using cluster 0')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stores last state for undo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    Sub_StoreUNDO
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data

% Increment Undo pointer - JCJ Oct 2007
if isfield(CBP_data,'UndoPtr') && ~isempty(CBP_data.UndoPtr)
    CBP_data.UndoPtr=CBP_data.UndoPtr+1;
else
    CBP_data.UndoPtr=1;
end
% Stuff UNDO
CBP_data.Undo(CBP_data.UndoPtr).a = CBP_data.a;
CBP_data.Undo(CBP_data.UndoPtr).h = CBP_data.h;
CBP_data.Undo(CBP_data.UndoPtr).iClust1 = CBP_data.iClust1;
CBP_data.Undo(CBP_data.UndoPtr).x1 = CBP_data.x1;
CBP_data.Undo(CBP_data.UndoPtr).y1 = CBP_data.y1;
CBP_data.Undo(CBP_data.UndoPtr).f1 = CBP_data.f1;
CBP_data.Undo(CBP_data.UndoPtr).iClust2 = CBP_data.iClust2;
CBP_data.Undo(CBP_data.UndoPtr).x2 = CBP_data.x2;
CBP_data.Undo(CBP_data.UndoPtr).y2 = CBP_data.y2;
CBP_data.Undo(CBP_data.UndoPtr).f2 = CBP_data.f2;
CBP_data.Undo(CBP_data.UndoPtr).keeps1 = CBP_data.keeps1;
CBP_data.Undo(CBP_data.UndoPtr).keeps2 = CBP_data.keeps2;
CBP_data.Undo(CBP_data.UndoPtr).Proj = CBP_data.Proj;

%if there are steps ahead of the new Undo step remove them -JCJ Oct 2007
if length(CBP_data.Undo)>CBP_data.UndoPtr
    CBP_data.Undo((CBP_data.UndoPtr+1):end)=[];
end

% Keep length of Undo cue limited to save memroy - JCJ Oct 2007
if (CBP_data.UndoPtr>CBP_data.MaxUndoNmbr)
    CBP_data.Undo(1)=[];
    CBP_data.UndoPtr=CBP_data.UndoPtr-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% moves to last state before last undo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sub_REDO()
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data

% Increment Undo pointer if further undo information exists and last non-undid state exists.
if ~isfield(CBP_data,'UndoPtr') || isempty(CBP_data.UndoPtr) || ~(CBP_data.UndoPtr<(length(CBP_data.Undo)-1))
    warndlg('No steps to redo.','CBP No Redo Available')
    return
else
    CBP_data.UndoPtr=CBP_data.UndoPtr+1;
end

% Use state before current undo (i.e. undo+one step).
CBP_Redo = CBP_data.Undo(CBP_data.UndoPtr+1);

CBP_data.a=CBP_Redo.a;
CBP_data.h=CBP_Redo.h;
CBP_data.iClust1=CBP_Redo.iClust1;
CBP_data.x1=CBP_Redo.x1;
CBP_data.y1=CBP_Redo.y1;
CBP_data.f1=CBP_Redo.f1;
CBP_data.iClust2=CBP_Redo.iClust2;
CBP_data.x2=CBP_Redo.x2;
CBP_data.y2=CBP_Redo.y2;
CBP_data.f2=CBP_Redo.f2;
CBP_data.keeps1=CBP_Redo.keeps1;
CBP_data.keeps2=CBP_Redo.keeps2;
CBP_data.Proj=CBP_Redo.Proj;

%Update display
set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',CBP_data.iClust1);
set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',CBP_data.iClust2 + 1);

clear CBP_Redo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Performs Undo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sub_RecallUndo()
% added restore of color information - JCJ Sept 2007
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data

% Decrement Undo pointer - JCJ Oct 2007
if ~isfield(CBP_data,'UndoPtr') || isempty(CBP_data.UndoPtr) || CBP_data.UndoPtr<1
    warndlg('No undo has been stored.','CBP No undo stored')
    return
end

% if stepping back for first time since an action, save current state ~ JCJ Oct 2007
if CBP_data.UndoPtr==length(CBP_data.Undo)
    OldPtr=CBP_data.UndoPtr;
    Sub_StoreUNDO;
    CBP_data.UndoPtr=OldPtr;
end

CBP_Undo = CBP_data.Undo(CBP_data.UndoPtr);
CBP_data.a=CBP_Undo.a;
CBP_data.h=CBP_Undo.h;
CBP_data.iClust1=CBP_Undo.iClust1;
CBP_data.x1=CBP_Undo.x1;
CBP_data.y1=CBP_Undo.y1;
CBP_data.f1=CBP_Undo.f1;
CBP_data.iClust2=CBP_Undo.iClust2;
CBP_data.x2=CBP_Undo.x2;
CBP_data.y2=CBP_Undo.y2;
CBP_data.f2=CBP_Undo.f2;
CBP_data.keeps1=CBP_Undo.keeps1;
CBP_data.keeps2=CBP_Undo.keeps2;
CBP_data.Proj=CBP_Undo.Proj;

%Update display
set(findobj('Tag', 'Cluster1pulldownCBP'),'Value',CBP_data.iClust1);
set(findobj('Tag', 'Cluster2pulldownCBP'),'Value',CBP_data.iClust2 + 1);

CBP_data.UndoPtr=CBP_data.UndoPtr-1; %JCJ Oct 2007

clear CBP_Undo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Clears Undo buffer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sub_ClearUNDO
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data

CBP_data.UndoPtr=0;
CBP_data.Undo=[];
warndlg('Cut On Best Projection: Undo buffer cleared.','CBP Undo Cleared')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Swaps info in Undo buffer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sub_SwapUNDO
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
global CBP_data

if isfield(CBP_data,'Undo');
    nUNDO=length(CBP_data.Undo);

    for iU=1:nUNDO
        Temp=CBP_data.Undo(iU);
        CBP_data.Undo(iU).h(1)=Temp.h(2);
        CBP_data.Undo(iU).h(2)=Temp.h(1);

        CBP_data.Undo(iU).iClust1=Temp.iClust2;
        CBP_data.Undo(iU).x1=Temp.x2;
        CBP_data.Undo(iU).y1=Temp.y2;
        CBP_data.Undo(iU).f1=Temp.f2;
        CBP_data.Undo(iU).keeps1=Temp.keeps2;

        CBP_data.Undo(iU).iClust2=Temp.iClust1;
        CBP_data.Undo(iU).x2=Temp.x1;
        CBP_data.Undo(iU).y2=Temp.y1;
        CBP_data.Undo(iU).f2=Temp.f1;
        CBP_data.Undo(iU).keeps2=Temp.keeps1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  Projection2dSuperposition(f1,f2,filename,ChannelValidity,iClust,iClust2,varargin)
% Projection2dSuperposition(f1,f2,filename,ChannelValidity,iClust);
%
% Projection2dSuperposition.m
%
% Finds a 2-D projection that is a superposition of the 3 best features to
% separate two clusters.
%
% INPUTS: f1,f2 -- index of points in the clusters 1 and 2
%         filename -- name of data file from which f was created
%         ChannelValidity -- 4 element vector specifying which channels to
%           use.
%         iClust  -- Cluster number to be cut on
%         iClust2 -- Cluster number to be separated
%
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
% JCJ 7 Dec 03 with strong input from ncst
% JCJ Oct 2004
% JCJ Feb 2008 -- included as sub-function of CutOnBestProjection
% STATUS: beta (in use)

% warning off MATLAB:divideByZero

if ~license('test', 'Statistics_Toolbox')
    warning('MCLUST:ToolboxUnavailable', 'Skipping cluster separation. Stats toolbox unavailable.');
    L_Extra = nan;
    L_Ratio = nan;
    IsolationDist = nan;
    Dists = [];
    df = nan;
    L_Intra = nan;
    C = nan;
    return
end

global CBP_data
global MClust_Colors MClust_Clusters 

ICASubsampThresh = 10000;% maximum number of spikes to include in ICA for one cluster -JCJ Sept 2007
MaxRecords = 400000;     % maximum number of records to load at one time
TTFileName = [];         % name of original tetrode file
n_records = [];          % number of records in TTFileName
FeaturesToGet = {};      % Features to use in calculating separation
FD = [];                 % Feature data matrix

if CBP_data.UseICA
    Nbest = 6;
else
    Nbest = 2;
end

extract_varargin;

if nargin < 3
    ChannelValidity = [1 1 1 1];
end

if nargin < 4
    iClust = ' ';
end

NewDir = 0;

[fpath fname fext] = fileparts(filename);

if isempty(TTFileName)
    TTFileName = [pwd filesep fname fext];
    if ~isempty(fpath)
        pushdir(fpath);
        TTFileName = [fpath filesep fname fext];
        NewDir = 1;
    end
end

if exist('FD','dir')
    FDdn = [pwd filesep 'FD'];
else
    FDdn = [pwd];
end

%%% Added to remove error when 'time' is selected as a feature
%%% JCJ Sept 2007
FeaturesToGet=FeaturesToGet(~strcmpi(FeaturesToGet,'time'));

nFeatures = length(FeaturesToGet);

FD1 = [];
FD2 = [];

pushdir(FDdn);
FDfiles = FindFiles([fname '_*.fd']);
popdir;

if ~isempty(FDfiles)
    temp = load(FDfiles{1},'-mat');
    if ~all(ChannelValidity == temp.ChannelValidity) & nargin < 3
        disp('Using Channel Validity from FD files')
        ChannelValidity = temp.ChannelValidity;
    end
end

% Find the size of the TT file
if isempty(n_records)
    n_records = MClust_CountSpikes(TTFileName);
end

% Length(f) needs to be longer than the number of columns in FD
if length(f1) < length(find(ChannelValidity))*length(FeaturesToGet) & length(f2) < length(find(ChannelValidity))*length(FeaturesToGet)
    disp('Not enough spikes to calculate mahalanobis distance')
    return
end

CalcFeatures = 0;
for iF = 1:nFeatures
    FDFileName = fullfile(FDdn, [fname '_' FeaturesToGet{iF} '.fd']);
    if ~exist(FDFileName)
        CalcFeatures = 1;
    end
end
if CalcFeatures
    Get_MClust_FD_Defaults
    Write_fd_file(FDdn, TTFileName, ...
        FeaturesToGet, ChannelValidity, record_block_size, template_matching, NormalizeFDYN)
end

if ~isempty(strmatch('keeps1',fields(CBP_data)))
    f1 = CBP_data.f1(find(CBP_data.keeps1));
    f2 = CBP_data.f2(find(CBP_data.keeps2));
    %     FD1 = FD1(find(CBP_data.keeps1),:);
    %     FD2 = FD2(find(CBP_data.keeps2),:);
end

FDnames={};
for iF = 1:nFeatures
    FDFileName = fullfile(FDdn, [fname '_' FeaturesToGet{iF} '.fd']);
    temp = load(FDFileName,'-mat');
    if size(temp.FeatureData,2) == size(ChannelValidity,2)
        FD1 = [FD1 temp.FeatureData(f1,find(ChannelValidity))];
        FD2 = [FD2 temp.FeatureData(f2,find(ChannelValidity))];
    else
        FD1 = [FD1 temp.FeatureData(f1,:)];
        FD2 = [FD2 temp.FeatureData(f2,:)];
    end
    FDnamesIdx=find(ChannelValidity); %find(strncmpi([FeaturesToGet{iF} ':'],MClust_FeatureNames,length(FeaturesToGet{iF})+1));
    for iC=1:length(FDnamesIdx)
        FDnames(end+1)=temp.FeatureNames(FDnamesIdx(iC)); %MClust_Temp_FeatureNames(FDnamesIdx(iC));
    end
end

if CBP_data.UseCQ & ~CBP_data.Proj.Save
    fdfn=findfiles('*.fd');
    pushdir(fileparts(fdfn{1}));

    [Best_x, Best_y] =FBP(FeaturesToGet,f1,f2,Nbest);

    for iB=1:Nbest
        ix(iB)=find(strcmpi(Best_x{iB},FDnames));
        iy(iB)=find(strcmpi(Best_y{iB},FDnames));
    end

    iD=unique([ix,iy]);

    popdir
else
    iD=1:length(FDnames);
end

if CBP_data.Proj.Save
    iD=CBP_data.Proj.iD;
else
    CBP_data.Proj.iD=iD;
end

FDnames=FDnames(iD);
FD1=FD1(:,iD);
FD2=FD2(:,iD);

lfd1 = size(FD1,1);
lfd2 = size(FD2,1);

%Added Sept 2007 JCJ (+ Changes to CutOnBestProjection and CBP_Callbacks)
if CBP_data.UseICA & ~CBP_data.Proj.Save

    % Sub-sample cluster 1 if needed (improves speed) Sept 2007 JCJ
    if lfd1>ICASubsampThresh
        idxFD1=randperm(lfd1);
        idxFD1=idxFD1(1:ICASubsampThresh);
    else
        idxFD1=1:lfd1;
    end

    % Sub-sample cluster 2 if needed (improves speed) Sept 2007 JCJ
    if lfd2>ICASubsampThresh
        idxFD2=randperm(lfd2);
        idxFD2=idxFD2(1:ICASubsampThresh);
    else
        idxFD2=1:lfd2;
    end
    % Put clusters together to define entire sample for ICA
    mixedsig=[FD1(idxFD1,:);FD2(idxFD2,:)]';


    figure(findobj('Name','Cut On Best Projection'))
    subplot('Position',CuttingPlotPosition)
    title('Please wait while idependent components are calculated...');
    pause(0.01)

    % Use ICA to pre-process feature data (tends to improve separation in
    % later stages).
    [A,W] = fastica (mixedsig,'numOfIC',length(iD),... %max(round(length(iD)/2),3),...
        'maxNumIterations',100,'verbose','off','stabilization','on');
    % now process full data set
    mixedsig=[FD1;FD2]';
    icasig = W*mixedsig;

    figure(findobj('Name','Cut On Best Projection'))
    subplot('Position',CuttingPlotPosition)
    title('ICA finished');

    % Separate samples corresponding to clusters
    FD1=icasig(:,1:lfd1)';
    FD2=icasig(:,(lfd1+1):end)';

    if (lfd1 ~= size(FD1,1))|(lfd2 ~= size(FD2,1))
        error('Incorrect sizes feature data after ICA')
    end

    CBP_data.ICA.W=W; % Save separating matrix for future use;

elseif CBP_data.Proj.Save
    mixedsig=[FD1;FD2]';
    % Use previously calculated ICA separating matrix reseparate feature data
    W= CBP_data.ICA.W;

    icasig = W*mixedsig;

    % Separate samples corresponding to clusters
    FD1=icasig(:,1:lfd1)';
    FD2=icasig(:,(lfd1+1):end)';

    if (lfd1 ~= size(FD1,1))|(lfd2 ~= size(FD2,1))
        error('Incorrect sizes feature data after ICA')
    end
else % Save identity matrix as separating matrix if ICA is not used
    %     (eliminates errors if CBP_data.Proj.Save is true);
    CBP_data.ICA.W=eye(length(iD));
end

if ~CBP_data.Proj.Save

    %first vector to define the plane is the vector between the two
    %cluster centers.
    mf1=mean(FD1);
    mf2=mean(FD2);
    V12x=(mf2-mf1)';

    %find the principal axes and eigenvalues for each cluster's
    %covariance matrix
    C1=cov(FD1);
    %C2=cov(FD2);

    %[EV1,ED1]=eigs(C1,length(C1));
    %[U,S,V,F]  =svds(C1,length(C1));
    %[V2,D2]=eigs(C2);
    [EV1,ED1,jnk,converge_flag]=svds(C1,length(C1));
    if converge_flag
        warning('Diagonalization did not converge, using reduced dimensionality')
        nD=length(C1);
        while converge_flag
            nD=nD-1;
            if nD>1
                [EV1,ED1,jnk,converge_flag]=svds(C1,nD);
            else
                error('Diagonalization not possible, use different feature space')
            end
        end
    end

    CBP_data.Proj.EV1=EV1;
    CBP_data.Proj.ED1=ED1;
    CBP_data.Proj.mf1=mf1;
else
    EV1=CBP_data.Proj.EV1;
    ED1=CBP_data.Proj.ED1;
    mf1=CBP_data.Proj.mf1;
    bestproj=CBP_data.Proj.iBP;
end


% Transform clusters into a coordinate system where cluster 1 is unit
% gaussian centered at origin and covariance matrix of cluster 1 is diagonal (D1)
c=((FD1-repmat(mf1,lfd1,1))*EV1)/sqrt(ED1);
d=((FD2-repmat(mf1,lfd2,1))*EV1)/sqrt(ED1);
% c=(EV1/sqrt(ED1))'*(FD1-repmat(mf1,length(FD1),1));
% d=(EV1/sqrt(ED1))'*(FD2-repmat(mf1,length(FD2),1));

if ~CBP_data.Proj.Save
    V12x=(EV1/sqrt(ED1))'*V12x; %vector from center of cluster1 to center of cluster2
    % V12x - (mean(d)-mean(c))' %test if transform worked. <-- it does to within round off ~1e-14
    % V12x = (median(d)-median(c))';

    C2=cov(d);  % covariance of cluster 2 in new coordinate system

    V12d2=inv(C2)*V12x.^2; % effective distance in each dimension from the cluster center

    % Use the 3 dimensions that provide the largest separation between
    % cluster centers given the distribution of cluster 2
    % [ys, yi]=sort(abs(V12x));
    [ys, yi] = sort(abs(V12d2));
    bestproj = yi(end-2:end);
    CBP_data.Proj.iBP=bestproj;
end

cb = c(:,bestproj);
db = d(:,bestproj);

if ~CBP_data.Proj.Save
    Vbx= V12x(bestproj);

    % Use the principal axis of Cluster 2 with the largest eigenvalue
    C2=cov(db);

    [temp1_val(1) temp1_idx(1)]=max(max(abs(C2)));
    %[temp1_val(1) temp1_idx(1)]=min(min(abs(C2)));
    Vbtemp=EV1(bestproj,temp1_idx(1));

    % Now calculate the plane of best projection for the 3 dimensions
    Vby=cross(Vbx,Vbtemp);

    CBP_data.Proj.Vby=Vby;
    CBP_data.Proj.Vbx=Vbx;
else
    Vby=CBP_data.Proj.Vby;
    Vbx=CBP_data.Proj.Vbx;
end

newFD1x=cb*Vbx;
newFD1y=cb*Vby;
newFD2x=db*Vbx;
newFD2y=db*Vby;

% figure
% hold on
% plot(newFD2x,newFD2y,'*','MarkerSize',MarkerSize,'Color',MClust_Colors(iClust2 + 1,:));
% plot(newFD1x,newFD1y,'*','MarkerSize',MarkerSize,'Color',MClust_Colors(iClust + 1,:));
idxS=ismember(f1,f2);
% pause

% figure
hold off

[Marker, MarkerSize] = GetMarker;

figure(findobj('Name','Cut On Best Projection'))
subplot('Position',CuttingPlotPosition)

if CBP_data.frontback
    CBP_data.h(1)=plot(newFD1x,newFD1y,Marker,'MarkerSize',MarkerSize,'Color',MClust_Colors(iClust + 1,:));
    hold on
    CBP_data.h(2)=plot(newFD2x,newFD2y,Marker,'MarkerSize',MarkerSize,'Color',MClust_Colors(iClust2 + 1,:));
else
    CBP_data.h(2)=plot(newFD2x,newFD2y,Marker,'MarkerSize',MarkerSize,'Color',MClust_Colors(iClust2 + 1,:));
    hold on
    CBP_data.h(1)=plot(newFD1x,newFD1y,Marker,'MarkerSize',MarkerSize,'Color',MClust_Colors(iClust + 1,:));
end
plot(newFD1x(idxS),newFD1y(idxS),'og','MarkerSize',MarkerSize+4);

if ~CBP_data.UseICA
    title(['Features used (best to worst): ' FDnames{bestproj(3)} ', ' FDnames{bestproj(2)} ', ' FDnames{bestproj(1)}])
end
a=axis;
% add a margin to make adding limits easier -- JCJ Sept 2007
da=[a(2)-a(1) a(4)-a(3)]/20;
a=a+[-da(1) da(1) -da(2) da(2)];
axis(a)
CBP_data.a=a;

if isempty(strmatch('keeps1',fields(CBP_data)))
    CBP_data.iClust1 = iClust;
    CBP_data.x1 = newFD1x;
    CBP_data.y1 = newFD1y;
    CBP_data.keeps1 = logical(repmat(1,length(newFD1x),1));
    CBP_data.f1 = f1;

    CBP_data.iClust2 = iClust2;
    CBP_data.x2 = newFD2x;
    CBP_data.y2 = newFD2y;
    CBP_data.keeps2 = logical(repmat(1,length(newFD2x),1));
    CBP_data.f2 = f2;
else
    f1_map = find(CBP_data.keeps1);
    CBP_data.x1(f1_map) = newFD1x;
    CBP_data.y1(f1_map) = newFD1y;

    f2_map = find(CBP_data.keeps2);
    CBP_data.x2(f2_map) = newFD2x;
    CBP_data.y2(f2_map) = newFD2y;
end

% Prepare UNDO (erase old UNDO data)
% CBP_data.Undo.x1 = CBP_data.x1;
% CBP_data.Undo.y1 = CBP_data.y1;
% CBP_data.Undo.x2 = CBP_data.x2;
% CBP_data.Undo.y2 = CBP_data.y2;
% CBP_data.Undo.keeps1 = CBP_data.keeps1;
% CBP_data.Undo.keeps2 = CBP_data.keeps2;

return

if NewDir
    popdir;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Marker, MarkerSize] = GetMarker
% [Marker, MarkerSize] = GetMarker
%
% Gets the current values [Marker, MarkerSize] of the
% Cluster Cutting Control Window plot settings and updates
% the global variables.
%
% JCJ Sept 2004


global MClust_ClusterCutWindow_Marker MClust_ClusterCutWindow_MarkerSize

figHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');
cboHandle = findobj('Parent',figHandle,'Tag', 'RedrawAxes');
figHandle = ParentFigureHandle(cboHandle);
markerHandle = findobj(figHandle, 'Tag', 'PlotMarker');
markerString = get(markerHandle, 'String');
markerValue = get(markerHandle, 'Value');
MClust_ClusterCutWindow_Marker = markerValue;
Marker = markerString{markerValue};

if nargout == 2
    markerSizeHandle = findobj(figHandle, 'Tag', 'PlotMarkerSize');
    markerSizeString = get(markerSizeHandle, 'String');
    markerSizeValue = get(markerSizeHandle, 'Value');
    MClust_ClusterCutWindow_MarkerSize = markerSizeValue;
    MarkerSize = str2num(markerSizeString{markerSizeValue});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

global MClust_Temp_FeatureNames MClust_Clusters  MClust_FDfn MClust_FeatureTimestamps

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

Curr_X = ' ';
Curr_Y = ' ';

DisplayProgress(0,nFeatures,'Title',['Finding Best Projections']);
for iX = 1:nFeatures
    FeatureToGet = MClust_Temp_FeatureNames{iX};
    FindColon = find(FeatureToGet == ':');
    if strmatch(upper(FeatureToGet(1:FindColon-1)),FBS_Features,'exact')
        if ~strcmpi(FeatureToGet(1:FindColon-1),Curr_X)
            temp_X = load(fullfile(fpath, [fname '_' FeatureToGet(1:FindColon-1) '.fd']),'-mat');
            Curr_X = FeatureToGet(1:FindColon-1);
        end
        FeatureIndex = strmatch(upper(FeatureToGet),upper(temp_X.FeatureNames));
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
                FeatureIndex = strmatch(upper(FeatureToGet),upper(temp_Y.FeatureNames)); %JCJ April 2008 -- convert to upper (case insensitive)
                FD_temp(:,2) = temp_Y.FeatureData(:,FeatureIndex);

                Dists = mahal(FD_temp(f_2,:),FD_temp(f_1,:));

                % Using a faster version of L
                [H bins] = hist(log10(Dists),500);
                Chi_CDF = 1 - chi2cdf(10.^bins,2);
                M_dists(iX,iY) = -sum(H.*Chi_CDF);
            end
        end;
    end
    DisplayProgress(iX,nFeatures);
end

Best_x = {}; Best_y = {};
for iP = 1:NBest
    [jnk yi] = max(max(M_dists));
    Best_y{iP} = MClust_Temp_FeatureNames{yi};
    [jnk xi] = max(M_dists(:,yi));
    Best_x{iP} = MClust_Temp_FeatureNames{xi};

    M_dists(xi,yi) = NaN;
end