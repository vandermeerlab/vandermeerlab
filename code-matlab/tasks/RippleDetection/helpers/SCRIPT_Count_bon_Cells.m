%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                        Count bon cells                              %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Goes through all bon data (rat Bond from Frank Lab) and counts the number
% of cells.
%
% bon data must be reorganized such that session data is grouped together
% in daily folders:
% D:\data\bon\bon-03
% D:\data\bon\bon-04
% D:\data\bon\bon-05
% D:\data\bon\bon-06
% D:\data\bon\bon-07
% D:\data\bon\bon-08
% D:\data\bon\bon-09
% D:\data\bon\bon-10
%
% In each folder, there must be a copy of boncellinfo.mat
%
% aacarey Dec 2017

area = 'CA1'; %'CA1','CA3', or '' for all areas

% Set current directory to the parent bon folder
cd('D:\data\bon') % bon datapath for VYSERITHUS

epoch = 7; % [] does all epochs

cfg = []; cfg.area = 'CA1';
cfg.verbose = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(epoch)
    epoch = 1:7;
end

cellCounts = nan(7,10); % store cells counts in an nEpoch x nSessions matrix

folderList = dir;
numFolders = length(folderList)-2;

bonDirectory = pwd;

disp(' ')
for iFolder = 3:length(folderList)
   cd([bonDirectory,'\',folderList(iFolder).name])
   recordingDay = iFolder; % this happens to be true for this rat. The dir list gives 2 empty results, and this rat is missing day 1 and day 2
   
   for iEpoch = 1:length(epoch)
       S = LoadFrankSpikes(cfg,epoch(iEpoch));
       numCells = length(S.t);
       cellCounts(epoch(iEpoch),recordingDay) = numCells;
       fprintf('Recording day %d, epoch %d: %d %s units\n',recordingDay,epoch(iEpoch),numCells,cfg.area)
   end
end

cd(bonDirectory)