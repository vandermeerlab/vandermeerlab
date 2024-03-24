%%
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

please = []; please.fc = ExpKeys.goodSWR(1); SWR_lfp_raw = LoadCSC(please);
please = []; please.fc = ExpKeys.goodTheta(1); theta_lfp_raw = LoadCSC(please);

pos = LoadPos([]); % also load position data (x,y) over time
%% get candidate events (SWR thresholding)
cfg = []; cfg.f = [140 220]; cfg.display_filter = 0;
SWR_lfp_filtered = FilterLFP(cfg, SWR_lfp_raw); theta_lfp_filtered = FilterLFP(cfg, theta_lfp_raw);
 
%% obtain differential power and z-score it
cfg = []; cfg.output = 'power';
SWR_power = LFPpower(cfg, SWR_lfp_filtered); control_power = LFPpower(cfg, theta_lfp_filtered);
diff_power_z = SWR_power; diff_power_z.data = diff_power_z.data - control_power.data; diff_power_z = zscore_tsd(diff_power_z);

%% set initial detection to be more lenient
cfg = []; cfg.method = 'raw'; cfg.threshold = 3; cfg.minlen = 0.025;
SWR_iv = TSDtoIV(cfg, diff_power_z);

% to each event, add a field with the max z-scored power (for later selection)
cfg = []; cfg.method = 'max'; cfg.label = 'maxSWRpower_z'; % what to call this in iv, i.e. usr.label
SWR_iv = AddTSDtoIV(cfg, SWR_iv, diff_power_z);
 
% select only those events of >5 z-scored power max
cfg = []; cfg.operation = '>'; cfg.threshold = 5;
SWR_iv = SelectIV(cfg,SWR_iv,'maxSWRpower_z');

%% set up place field templates for left and right trials
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
   
    this_coord.coord = expCond(iCond).coord; this_coord.units = 'au'; this_coord.standardized = 0;
    expCond(iCond).linpos = LinearizePos([],pos,this_coord);
   
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

%% get tuning curves and associated template (order by field location)
for iCond = 1:nCond
    cfg_tc = [];
    expCond(iCond).tc = TuningCurves(cfg_tc,expCond(iCond).S,expCond(iCond).linpos); 
    expCond(iCond).fields = DetectPlaceCells1D([],expCond(iCond).tc.tc);
end

figure; 

subplot(221)
imagesc(expCond(1).tc.tc); xlabel('linearized position (bin)'); ylabel('cell ID'); title('left');
subplot(223)
left_template = expCond(1).fields.template_idx;
imagesc(expCond(1).tc.tc(left_template, :)); xlabel('linearized position (bin)'); ylabel('template cell ID');

subplot(222)
imagesc(expCond(2).tc.tc); xlabel('linearized position (bin)'); ylabel('cell ID'); title('right');
subplot(224)
right_template = expCond(2).fields.template_idx;
imagesc(expCond(2).tc.tc(right_template, :)); xlabel('linearized position (bin)'); ylabel('template cell ID');


%% method for replay scoring: correlate template order with spike time order
addpath(genpath(cat(2,SET_GitHub_root,'\vandermeerlab\code-matlab\tasks\Replay_Analysis'))); % clone vandermeerlab repo at https://github.com/vandermeerlab/vandermeerlab

% note: could improve candidate selection here (e.g. animal stationary, not
% in theta, minimum number of cells active, etc.)

S_left = SelectTS([],S,expCond(1).fields.template_idx); 
cfg = [];
cfg.nShuffles = 100;
cfg.seqMethod = 'exact';
score_left = scoreCandSeq(cfg, SWR_iv, S_left); %% TODO -- check if this works with improved function, YT's CorrScoreWin3()?

% add scores to usr field of candidate ivs. N.B.- this can be done in scoreCandSeq() but
% for now we have a "scores" output which we can append back onto the candidate ivs.
datas = struct2cell(score_left);
label = fieldnames(score_left);
for f = 1:length(label)
    SWR_iv.usr(f).data = datas{f}';
    SWR_iv.usr(f).label = label{f};
end

%% display results
cfg = [];
cfg.lfp = SWR_lfp_raw; cfg.lfpColor = 'k';
cfg.evt = SWR_iv;
MultiRaster(cfg,S_left);
ylim([-25 45])