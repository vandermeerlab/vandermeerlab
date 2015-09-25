%% preamble
clear all; pack

%MASTER_root = 'D:\My_Documents\Dropbox\projects\Alyssa'; % bergkamp
MASTER_root = 'C:\Users\mvdm\Dropbox\projects\Alyssa'; % isidro
cd(MASTER_root);

MASTER_path;

%% set up fds
fd_cfg = []; fd_cfg.requireCandidates = 0;
fd = getTmazeDataPath(fd_cfg);

cd(fd{end});

%% load data (quickload?, FileLoader{'spikes','pos','ExpKeys','Metadata'))
S = LoadSpikes([]);
pos = LoadPos([]);
LoadExpKeys;
LoadMetadata;

%% specify experimental conditions
clear expCond;

expCond(1).label = 'left';
expCond(2).label = 'right';

expCond(1).t = metadata.taskvars.trial_iv_L;
expCond(2).t = metadata.taskvars.trial_iv_R;

expCond(1).coord = metadata.coord.coordL; % note, pizels to match loadpos()
expCond(2).coord = metadata.coord.coordR; % note, pizels to match loadpos()

%% linearize paths
nCond = length(expCond);
for iCond = 1:nCond
   
    cfg_linpos = []; 
    cfg_linpos.Coord = expCond(iCond).coord;
    expCond(iCond).linpos = LinearizePos(cfg_linpos,pos);
   
end

%% get speed
spd = getLinSpd([],pos);

cfg_spd = []; 
cfg_spd.method = 'raw'; 
cfg_spd.threshold = 10; 
run_iv = TSDtoIV(cfg_spd,spd);

%% set up data to use for TC
for iCond = 1:nCond
   
    expCond(iCond).S = restrict(S,run_iv);
    expCond(iCond).S = restrict(S,expCond(iCond).t);

    % modify this using structfun (needs to work on specific fields)
    expCond(iCond).linpos = restrict(expCond(iCond).linpos,run_iv);
    expCond(iCond).linpos = restrict(expCond(iCond).linpos,expCond(iCond).t);
    
end

%% get tuning curves
for iCond = 1:nCond
    
    cfg_tc = [];
    expCond(iCond).tc = TuningCurves(cfg_tc,expCond(iCond).S,expCond(iCond).linpos); 
    
end

%% find some fields
for iCond = 1:nCond
    
    expCond(iCond).fields = DetectPlaceCells1D([],expCond(iCond).tc.tc);
    
end
    
%% load CSC
cfg = [];
cfg.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(cfg);

%% plot something
fh = figure('KeyPressFcn',@navigate);

for iCond = 1:nCond
ax(iCond) = subplot(2,1,iCond);
cfg_mr = []; cfg_mr.openNewFig = 0; cfg_mr.lfp = csc;
S_temp = S; S_temp.t = S_temp.t(expCond(iCond).fields.template_idx);
MultiRaster(cfg_mr,S_temp);
end
linkaxes(ax,'x')