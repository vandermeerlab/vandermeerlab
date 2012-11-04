function [redraw, rekey, undoable] = WriteDotCutFile

% INPUTS
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
%
% ADR 2010
%
% Status: PROMOTED (Release version)
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

redraw = false; rekey = false; undoable = false;

global MClust_Clusters MClust_FeatureTimestamps MClust_TTfn  MClust_TTdn

dotCutFilename = fullfile(MClust_TTdn, [MClust_TTfn '.cut']);
dotCutFilehandle = fopen(dotCutFilename, 'wt');

nS = length(MClust_FeatureTimestamps);
nClust = length(MClust_Clusters);

spikeClusterMatrix = false(nClust, nS);
for iC = 1:nClust
    [fI] = FindInCluster(MClust_Clusters{iC});
    spikeClusterMatrix(iC, fI) = true;
end
for iS = 1:nS
    fprintf(dotCutFilehandle,'%d: ', iS);
    for iC = 1:nClust
        if spikeClusterMatrix(iC, iS)
            fprintf(dotCutFilehandle,'%d ', iC);
        end
    end
    fprintf(dotCutFilehandle,'\n');
end

fclose(dotCutFilehandle);
