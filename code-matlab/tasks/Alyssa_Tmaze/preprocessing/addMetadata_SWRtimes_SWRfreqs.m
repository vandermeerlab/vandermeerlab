% Get SWRfreaks (Fourier coefficients) and SWRtimes for T-maze metadata

clear

% First get a set of SWRfreqs from a past session to help find good SWRs more
% quickly:
originalFolder = pwd;
 cd('D:\data\R050\R050-2014-04-02')

 LoadMetadata

SWRfreqs_temp = metadata.SWRfreqs;

clear metadata

cd(originalFolder)

cfg = [];
cfg.useClustersFile = 0;
cfg.load_questionable_cells = 1;
S = LoadSpikes(cfg);

LoadExpKeys
cfg.fc = ExpKeys.goodSWR(1);
%cfg.fc = {'R044-2013-12-22-CSC02a.ncs'};
csc = LoadCSC(cfg);

[SWR,~,~] = amSWR([],SWRfreqs_temp,csc);
[MUA,~,~] = amMUA([],S,csc.tvec);
cfg_can.threshold = 13; % set threshold pretty high so don't have to move through false positives
evt = precand(cfg_can,csc.tvec,SWR,MUA,S);
%evt = detector([],SWRfreqs_temp,csc,S);

%% manually identify SWRs 

cfg_duck.num = 50;
cfg_duck.evt = evt;
cfg_duck.lfpHeight = 30;

SWRtimes = ducktrap(cfg_duck,S,csc);

%test.tcent = ducktrap_backup;
%test.label = 'blah';

%% get SWR Fourier coefficients
cfg.showfig = 1;
cfg.weightby = 'amplitude';
SWRfreqs = SWRfreak(cfg,SWRtimes,csc);

%% Save fields in metadata

% WARNING: running this section overwrites existing SWRtimes and SWRfreqs fields (if they exist already)!

% first check if metadata exists yet

loaded = LoadMetadata2;

if ~loaded % then it doesn't exist yet, so make a metadata struct
    metadata.SWRtimes = SWRtimes;
    metadata.SWRfreqs = SWRfreqs;
else % it does, so add a new field
    metadata.SWRtimes = SWRtimes;
    metadata.SWRfreqs = SWRfreqs;
end

[~,name,~] = fileparts(pwd); % pwd is your current folder, we just want its namepart

% now save
savename = strcat(name,'-metadata.mat'); % use the folder's name, but concatenate it with '-metadata'
save(savename,'metadata'); % this saves the specified variables under the given [save]name