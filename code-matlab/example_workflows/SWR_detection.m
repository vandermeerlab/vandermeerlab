%% setup
SET_GitHub_root = 'C:\Users\mattm\Documents\GitHub'; % replace this with the location of your local GitHub root folder
%SET_GitHub_root = 'C:\Users\mvdm\Documents\GitHub';
SET_data_fd = 'C:\data\R064\R064-2015-04-20'; % replace this with the location of your local data folder, download from https://drive.google.com/file/d/19lO0KQQjK9NQjr3xOz38bHwaG2yctr_q/view?usp=drive_link

restoredefaultpath;
addpath(genpath(cat(2,SET_GitHub_root,'\vandermeerlab\code-matlab\shared'))); % clone vandermeerlab repo at https://github.com/vandermeerlab/vandermeerlab

cfg_master = []; % holds config options for later in this script

%% load some data
cd(SET_data_fd);
 
please = []; please.load_questionable_cells = 1;
S = LoadSpikes(please); % `please` variable overrides default LoadSpikes() options to load all cells
 
LoadExpKeys; % annotation file containing some basic information about this data session
LoadMetadata; % loads experimenter-annotated file associated with each data session

please = []; please.fc = ExpKeys.goodSWR(1);
SWR_lfp_raw = LoadCSC(please);

please = []; please.fc = ExpKeys.goodTheta(1);
theta_lfp_raw = LoadCSC(please);

%% observe some candidate SWRs
cfg_mr = [];
cfg_mr.lfp(1) = SWR_lfp_raw; cfg_mr.lfp(2) = theta_lfp_raw;
MultiRaster(cfg_mr, S);

xlim([7005.5 7008.5]); % example of some SWRs, for session C:\data\R064\R064-2015-04-20

%% get ripple-band power
% filter in SWR band for rodents, 140-220 Hz
cfg = []; cfg.f = [140 220]; cfg.display_filter = 0;
SWR_lfp_filtered = FilterLFP(cfg, SWR_lfp_raw);
theta_lfp_filtered = FilterLFP(cfg, theta_lfp_raw);
 
%% obtain power and z-score it
cfg = []; cfg.output = 'power';
SWR_power = LFPpower(cfg, SWR_lfp_filtered); SWR_power_z = zscore_tsd(SWR_power);
control_power = LFPpower(cfg, theta_lfp_filtered); control_power_z = zscore_tsd(control_power);

% also compute a z-scored "power difference" of control channel ripple
% power subtracted from SWR channel ripple power -- this is a good way to
% reduce confounds from movement/other artifacts that are common across
% channels
diff_power_z = SWR_power; diff_power_z.data = diff_power_z.data - control_power.data; diff_power_z = zscore_tsd(diff_power_z);

%% plot z-scored ripple power alongside 
cfg_mr = []; cfg_mr.lfpMax = Inf; cfg_mr.lfpHeight = 15; cfg_mr.lfpSpacing = 5;
cfg_mr.lfp(1) = SWR_lfp_raw; cfg_mr.lfp(2) = SWR_power_z; cfg_mr.lfp(3) = diff_power_z;
MultiRaster(cfg_mr, S);

xlim([7005.5 7008.5]); % example of some SWRs, for session C:\data\R064\R064-2015-04-20

%% now threshold
cfg = []; cfg.method = 'raw'; cfg.threshold = 5;
SWR_iv = TSDtoIV(cfg, diff_power_z);

cfg_mr = []; cfg_mr.lfpMax = Inf; cfg_mr.lfpHeight = 15; cfg_mr.lfpSpacing = 5; cfg_mr.lfpColor = 'k';
cfg_mr.lfp(2) = SWR_lfp_raw; cfg_mr.lfp(1) = diff_power_z; cfg_mr.evt = SWR_iv;
MultiRaster(cfg_mr, S);

xlim([7005.5 7008.5]);

% notice we only detected one of the SWRs -- how to improve?
%% set initial detection to be more lenient
cfg = []; cfg.method = 'raw'; cfg.threshold = 3; cfg.minlen = 0.025;
SWR_iv = TSDtoIV(cfg, diff_power_z);

% to each event, add a field with the max z-scored power (for later selection)
cfg = [];
cfg.method = 'max'; % 'min', 'mean'
cfg.label = 'maxSWRpower_z'; % what to call this in iv, i.e. usr.label
SWR_iv = AddTSDtoIV(cfg, SWR_iv, diff_power_z);
 
%% select only those events of >5 z-scored power max
cfg = [];
cfg.operation = '>';
cfg.threshold = 5;
SWR_iv = SelectIV(cfg,SWR_iv,'maxSWRpower_z');

%% plot
cfg_mr.evt = SWR_iv;
MultiRaster(cfg_mr, S);

xlim([7005.5 7008.5]); 
% better, we're detecting both. Now look around and see how this detector
% performs for other candidate events -- what do you observe? can you
% tweak parameters more?

%% plot events only (slow)
cfg = []; cfg.display = 'iv'; cfg.mode = 'center'; cfg.fgcol = 'k';
PlotTSDfromIV(cfg, SWR_iv, SWR_lfp_raw);

%hold on; % optional: overlay detected event in red
%cfg = []; cfg.display = 'iv'; cfg.fgcol = 'r';
%PlotTSDfromIV(cfg, SWR_iv, SWR_lfp_raw);

%% different approach: curate test set, then train classifier
addpath(genpath(cat(2,SET_GitHub_root,'\vandermeerlab\code-matlab\tasks\Alyssa_Tmaze')));

please = [];
please.evt = SWR_iv;
please.mode = 'unfixed';
ducktrap(please, S, SWR_lfp_raw)
ylim([-25 116])

%% once training set is saved, run detector based on similarity to training examples
load(FindFile('*manualIV.mat'));
ncfs = SWRfreak([], evt, SWR_lfp_raw);

SWR_score = amSWR([], ncfs, SWR_lfp_raw); % might take ~1 minute
SWR_score_z = zscore_tsd(SWR_score);

%% plot to compare different SWR scoring methods
cfg_mr = []; cfg_mr.lfpMax = Inf; cfg_mr.lfpHeight = 15; cfg_mr.lfpSpacing = 5; %cfg_mr.lfpColor = 'kr';
cfg_mr.lfp(1) = SWR_score_z; 
diff_power_z.parameters = []; cfg_mr.lfp(2) = diff_power_z;
MultiRaster(cfg_mr, S);

xlim([7005.5 7008.5]);
