function [redraw, rekey, undoable]=SortClustersByColor

% [redraw, rekey, undoable] = SortClustersByColor
%
% SortClustersByColor  maps colors onto a circular color
% wheel, then sorts them by phase and length of the radial
% vector mapping.
%
% INPUTS
%   
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/GeneralizedCutterOptions folder

% JCJ Sept 2007
% JCJ April 2008 Updated for MClust3.5 compatibility

% JCJ April 2008 Added for MClust3.5 compatibility
redraw   = true;
rekey    = true;
undoable = true;

global MClust_Clusters MClust_Colors MClust_Hide
% JCJ April 2008 Removed MClust_ClusterFileNames for MClust3.5 compatibility

nClust  =length(MClust_Clusters);

[P, A]=Sub_Color2Phase(MClust_Colors(2:nClust+1,:))

[uColors, iColor, jColor] = unique([P A],'rows');

[sColor, Cidx]=sort(jColor);

% JCJ April 2008 Removed for MClust3.5 compatibility
% GeneralizedCutterStoreUndo('Sort By Color'); 

% MClust_Colors  =uColors(jColor(Cidx),:);
MClust_Colors(2:(nClust+1),:)  =MClust_Colors(Cidx+1,:);
MClust_Hide(2:nClust+1)    =MClust_Hide(Cidx+1);
MClust_Clusters=MClust_Clusters(Cidx);

% JCJ April 2008 Removed for MClust3.5 compatibility
% MClust_ClusterFileNames=MClust_ClusterFileNames(Cidx); 
% 
% ClusterCutWindow = findobj('Type','figure','Tag', 'ClusterCutWindow');
% if ~isempty(ClusterCutWindow)
%     close(ClusterCutWindow);
% end
% 
% CHDrawingAxisWindow = findobj('Type','figure','Tag', 'CHDrawingAxisWindow');
% if ~isempty(CHDrawingAxisWindow)
%     close(CHDrawingAxisWindow);
% end
% 
% MClustCallbacks CutPreClusters 

warndlg(['Clusters sorted by color.'], 'Color sort successful. Sort is Undo-able');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% This function maps colors onto a circular color wheel, then finds the
%%%%% phase in radians blue is equivalent to straight down (i.e. 270deg).
%%%%% C is a color [r g b]
%%%%% P is a phase in radians
%%%%% A is the amplitude
%%%%% JCJ Sept 2007;
function [P, A]=Sub_Color2Phase(C)
x=C*[2/sqrt(3); -2/sqrt(3); 0];
y=C*[1/sqrt(3); 1/sqrt(3); -1];
P=atan2(y,x);
P(P<0)=pi-P(P<0); % unwrap the vector
A=sqrt(x.^2+y.^2);