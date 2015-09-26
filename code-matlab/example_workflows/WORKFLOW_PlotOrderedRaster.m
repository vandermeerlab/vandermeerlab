%% preamble -- see lab wiki for explanation at:
% http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week1
clear all; pack

MASTER_root = 'D:\My_Documents\Dropbox\projects\Alyssa'; % replace this with home folder of your project
cd(MASTER_root);

MASTER_path; % reset and then set up path for this project

%% find data folders
fd_cfg = []; fd_cfg.requireCandidates = 0; % don't need previously saved sharp wave-ripple (SWR) candidates here
fd = getTmazeDataPath(fd_cfg);

cd(fd{end});

%% load data (alternative: FileLoader{'spikes','pos','ExpKeys','Metadata')), see lab wiki for loader info at:
% http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week2
S = LoadSpikes([]);
pos = LoadPos([]); % note, this is raw position data read from the .nvt file, units are "camera pixels"
LoadExpKeys; % see https://github.com/mvdm/vandermeerlab/blob/master/doc/HOWTO_ExpKeys_Metadata.md
LoadMetadata;

%% set up data structs for 2 experimental conditions -- see lab wiki for this task at:
% http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:task:motivationalt
clear expCond;

expCond(1).label = 'left'; % this is a T-maze, we are interested in 'left' and 'right' trials
expCond(2).label = 'right'; % these are just labels we can make up here to keep track of which condition means what

expCond(1).t = metadata.taskvars.trial_iv_L; % previously stored trial start and end times for left trials
expCond(2).t = metadata.taskvars.trial_iv_R; 

expCond(1).coord = metadata.coord.coordL; % previously user input idealized linear track
expCond(2).coord = metadata.coord.coordR; % note, this is in units of "camera pixels", not cm

expCond(1).S = S;
expCond(2).S = S;

%% linearize paths (snap x,y position samples to nearest point on experimenter-drawn idealized track)
nCond = length(expCond);
for iCond = 1:nCond
   
    cfg_linpos = []; cfg_linpos.Coord = expCond(iCond).coord;
    expCond(iCond).linpos = LinearizePos(cfg_linpos,pos);
   
end

%% find intervals where rat is running
spd = getLinSpd([],pos); % get speed (in "camera pixels per second")

cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 10; 
run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above 10 pix/s

%% restrict (linearized) position data and spike data to desired intervals
for iCond = 1:nCond
   
    fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
end

%% get tuning curves, see lab wiki at:
% http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week12
for iCond = 1:nCond
    
    cfg_tc = [];
    expCond(iCond).tc = TuningCurves(cfg_tc,expCond(iCond).S,expCond(iCond).linpos); 
    
end
 
%% detect place fields
for iCond = 1:nCond
    
    expCond(iCond).fields = DetectPlaceCells1D([],expCond(iCond).tc.tc);
    
end
    
%% load CSC with good SWRs
cfg = []; cfg.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(cfg);

%% make rasterplot place cells for left and right separately, ordered by place field location
fh = figure('KeyPressFcn',@navigate);

for iCond = 1:nCond
    ax(iCond) = subplot(2,1,iCond);
    
    S_temp = S; S_temp.t = S_temp.t(expCond(iCond).fields.template_idx); % template_idx contains ordered place cells
    
    cfg_mr = []; cfg_mr.openNewFig = 0; cfg_mr.lfp = csc;
    MultiRaster(cfg_mr,S_temp); % see http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week3short
end
linkaxes(ax,'x')