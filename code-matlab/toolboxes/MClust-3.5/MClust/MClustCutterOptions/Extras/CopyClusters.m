function [redraw, rekey, undoable]=Copy
% Copy
%  Copies clusters specified by user in a popup window.
%
% INPUTS: None (user interface)
%
% OUTPUTS: None (modifies MClust global variables)
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/GeneralizedCutterOptions folder

% JCJ Sept 2007
% JCJ April 2008 -- Updated for MClust 3.5
global MClust_Clusters MClust_ClusterFileNames MClust_Colors MClust_Hide
redraw   = true;
rekey    = true;
undoable = true;

figHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');

prompt={'Cluster(s) to Copy (comma or space separated list):','Copy Hides:','Copy Shows:'};
def={'','no','no'};
dlgTitle='Copy Cluster(s)';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
if ~isempty(answer)
    iClust1=str2num(answer{1});

    iClust2=[];
    if strncmpi(answer{2},'yes',1) | str2num(answer{2})==1
        nClust=length(MClust_Clusters);
        iClust2=find(MClust_Hide(2:(nClust+1)));
    end

    iClust3=[];
    if strncmpi(answer{3},'yes',1) | str2num(answer{3})==1
        nClust=length(MClust_Clusters);
        iClust3=find(~MClust_Hide(2:(nClust+1)));
    end

    iClust2Copy=unique([shiftdim(iClust1);iClust2;iClust3]);

    if ~isempty(iClust2Copy)
        nClust2Copy=length(iClust2Copy);

        for iD=1:nClust2Copy
            iClust=iClust2Copy(iD);
            
            MClust_Clusters{end+1} = MClust_Clusters{iClust};
            MClust_Hide(end+1) = 0;
            MClust_Colors(length(MClust_Clusters)+1,:) = MClust_Colors(iClust+1,:);
           
        end
        warndlg(['Cluster(s) ' num2str(iClust2Copy') ' copied.'], 'Copy successful. Undo Saved');
    else
        warndlg('No clusters copied.', 'Copy unsuccessful.');
    end
end
