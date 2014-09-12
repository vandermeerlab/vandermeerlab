function [redraw, rekey, undoable]=EvalOverlap(iClust)
redraw    = false;
rekey     = false;
undoable  = false;

global MClust_Clusters MClust_FeatureTimestamps
prompt={'Cluster to Compare: (space separated list)'};
def={''};
dlgTitle=['Eval Overlap with Cluster ' num2str(iClust)];
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
temp=str2num(answer{1});
nClust=length(temp)+1;
iClust(2:nClust)=temp;
    
    [nS, nD] = size(MClust_FeatureTimestamps);
    %nClust = length(MClust_Clusters);
    
    % ncst added 16 may 02
    if nClust == 1
        errordlg('Only one cluster exists.', 'MClust error', 'modal');
        return
    end
    
    nToDo = nClust * (nClust-1)/2;
    iDone = 0;
    overlap = zeros(nClust);
    for iInd = 1:nClust
        iC=iClust(iInd);
        iSpikes = zeros(nS, 1);
        [fI MClust_Clusters{iC}] = FindInCluster(MClust_Clusters{iC});
        iSpikes(fI) = 1;
        for jInd = (iInd+1):nClust
            jC=iClust(jInd);
            iDone = iDone +1;
            DisplayProgress(iDone, nToDo, 'Title', 'Evaluating overlap');
            jSpikes = zeros(nS, 1);
            [fJ MClust_Clusters{jC}] = FindInCluster(MClust_Clusters{jC});
            jSpikes(fJ) = 1;
            overlap(iInd,jInd) = sum(iSpikes & jSpikes);
            overlap(jInd,iInd) = overlap(iInd,jInd);
        end
    end
    overlap = [([0 iClust])', [iClust; overlap]]
    msgbox(num2str(overlap), 'Overlap');
    
