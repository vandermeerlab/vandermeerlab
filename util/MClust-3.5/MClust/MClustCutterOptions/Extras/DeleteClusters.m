function  [redraw, rekey, undoable]=Delete

% Delete
%  Deletes clusters specified by user in a popup window.
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


prompt={'Cluster(s) to Delete (comma or space separated list):','Delete Hides:','Delete Shows:'};
def={'','no','no'};
dlgTitle='Delete Cluster(s)';
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

    iClust2Delete=unique([shiftdim(iClust1);iClust2;iClust3]);

    if ~isempty(iClust2Delete)
        %         GeneralizedCutterStoreUndo('Delete Cluster');
        nClust2Delete=length(iClust2Delete);

        for iD=1:nClust2Delete
            iClust=iClust2Delete(iD);
            if isa(MClust_Clusters{iClust},'mccluster')
                [GetInfoMSG] = GetInfo(MClust_Clusters{iClust});
                if ~strcmpi(GetInfoMSG{strmatch('Limits',GetInfoMSG)+1},'none')
                    %                     for iX = 1:length(ClusterLimits(:,1))
                    MClust_Clusters{iClust} = CTO_11_Delete_All_Limits(MClust_Clusters{iClust});
                    %                     end
                end
                MClust_Clusters{iClust} =  CTO_13_Delete_All_Points(MClust_Clusters{iClust});

            end
        end
        %         GeneralizedCutterCallbacks('RedrawAxes')
        warndlg(['Cluster(s) ' num2str(iClust2Delete') ' deleted.'], 'Delete successful. Undo Saved');
    else
        warndlg('No clusters deleted.', 'Delete unsuccessful.');
    end
end
