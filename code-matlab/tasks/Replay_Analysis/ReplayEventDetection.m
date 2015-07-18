%% Load data
clear all; clc;
cd('D:\data\inProcess\R042-2013-08-18_recording'); %your data path
profile on;

% load a CSC 
cfg = []; cfg.fc = {'R042-2013-08-18_recording-CSC11b.Ncs'}; %whatever you call it
lfp = LoadCSC(cfg);

% Load spikes
cfg = [];
cfg.useClustersFile = 0;
S_orig = LoadSpikes(cfg);
t = [S_orig.cfg.ExpKeys.TimeOnTrack S_orig.cfg.ExpKeys.TimeOffTrack];
S_res = restrict2(S_orig,t(1),t(2)); %take out pre/post record

% Load Position data
cfg = [];
cfg.realTrackDims = [167 185];
pos = LoadPos(cfg);

% detect trial types
cfg = [];
cfg.eventList = {'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0020).',...
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0080).',...
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).',...
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0004).',...
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0040).'};
cfg.eventLabel = {'center_pb','left_pb','right_pb','food_dispensed','water_dispensed'};
evt = LoadEvents(cfg);

load('times.mat');
evt.run_start = run_start;
evt.run_end = run_end;
evt.ped_start = ped_start*10^-6;
evt.ped_end = ped_end*10^-6;
clear nBins ped_end ped_start prerecord_end prerecord_start run_end run_start t_end t_start

% separate left/right trials 
for iRun = 1:length(evt.run_start)
    temp = getd(evt,'food_dispensed');
    if ~isempty(temp(temp > evt.run_start(iRun) & temp < evt.run_end(iRun)))
        evt.trial_ID(iRun) = 1; % food/left
    else
        evt.trial_ID(iRun) = 2; % water/right
    end
end

left_trials = find(evt.trial_ID == 1);
evtL = evt;
evtL.trial_ID = evtL.trial_ID(left_trials);
evtL.run_start = evtL.run_start(left_trials);
evtL.run_end = evtL.run_end(left_trials);

right_trials = find(evt.trial_ID == 2);
evtR = evt;
evtR.trial_ID = evtR.trial_ID(right_trials);
evtR.run_start = evtR.run_start(right_trials);
evtR.run_end = evtR.run_end(right_trials);
clear temp iRun left_trials right_trials

% Separate left/right spikes & positions
S_left = restrict(S_res,evtL.run_start,evtL.run_end);
posL = restrict(pos,evtL.run_start,evtL.run_end);
S_right = restrict(S_res,evtR.run_start,evtR.run_end);
posR = restrict(pos,evtR.run_start,evtR.run_end);

% Linearize runs and bin positions
load('D:\data\pfmodels\2014-12-21\coord.mat'); % Use MakeCoord if you don't already have this
% CoordL = MakeCoord(getd(posL,'x'),getd(posL,'y'));
% CoordR = MakeCoord(getd(posR,'x'),getd(posR,'y'));

run_dist = 274; % distance travelled on a single run of the track in cm (T-maze)
binSize = 2; % in cm (2 for D&G, 1 for Davidson et al)
nBins = run_dist/binSize;
CoordLrs(1,:) = interp1(1:size(CoordL,2),CoordL(1,:),linspace(1,size(CoordL,2),nBins),'linear');
CoordLrs(2,:) = interp1(1:size(CoordL,2),CoordL(2,:),linspace(1,size(CoordL,2),nBins),'linear');
CoordRrs(1,:) = interp1(1:size(CoordR,2),CoordR(1,:),linspace(1,size(CoordR,2),nBins),'linear');
CoordRrs(2,:) = interp1(1:size(CoordR,2),CoordR(2,:),linspace(1,size(CoordR,2),nBins),'linear');

cfg = [];
cfg.Coord = CoordLrs;
posL_binned = LinearizePos(cfg,posL);
cfg = [];
cfg.Coord = CoordRrs;
posR_binned = LinearizePos(cfg,posR);

%store left/right linearized data in a single struct (cleaner workspace!)

data(1).trial_type = 'left';
data(1).evt = evtL;
data(1).Coord = CoordLrs;
data(1).pos = posL_binned;
data(1).S = S_left;

data(2).trial_type = 'right';
data(2).evt = evtR;
data(2).Coord = CoordRrs;
data(2).pos = posR_binned;
data(2).S = S_right;

clear evtL evtR CoordL CoordLrs CoordR CoordRrs S_left S_right posL posL_binned posR posR_binned nBins trial_ID t
 
%% Make tuning curves for left or right trajectorie 
clear TC %in case you are re-running
cfg = [];
cfg.binSize = binSize;
iT = 1; %1 for left / 2 for right

TC = MakeTC2(cfg,data(iT).S,data(iT).pos);

%% Plot tuning curves - visualization like Dragoi & Tonegawa (2011)
figure;
tc_temp = TC.tc(:,TC.template_idx)';
for iC = 1:size(tc_temp,1)       
    subaxis(size(tc_temp,1),1,iC,'SpacingVert',0);
    area(tc_temp(iC,:)); hold on;
    ylims = get(gca,'Ylim');
    plot([TC.peak_idx(iC) TC.peak_idx(iC)],[ylims(1) ylims(2)],'r:');
    set(gca, 'XTick', [],'YTick',[]); 
    yl = ylabel(iC); 
    set(yl,'Rotation',0,'Fontsize',8);
    if iC == 1; title(data(iT).trial_type); end
end
x = size(tc_temp,2);
set(gca,'XTick',[1 ceil(x/2) x],'XTicklabel',[1 137 274],'Ticklength', [0 0],'Fontsize',8);

clear iC tc_temp trial_ID x yl ylims

%% Designate place cells
field_order = TC.template_idx;
S_pc = S_orig;
S_pc.t = S_pc.t(field_order); %arrange spike train by place cell ordering
S_pc.label = S_pc.label(field_order);

%% Exract SWR events
cfg = [];
cfg.threshold = 3;
cfg.max_thr = 3;
spwr_iv = getSWR(cfg,S_pc,lfp);

%% plot SWR events
cfg = [];
cfg.lfp = lfp;
cfg.evt = spwr_iv;
cfg.windowSize = 0.5; %in seconds
MultiRaster(cfg,S_pc); %NOTE: MultiRaster has the navigate function inside!!!

%% Create candidate events
% SWR events are passed into getCandSeq as a cfg parameter. getSWR can easily be called
% inside getCandSeq, but this way let's us examine the SWR separately.

cfg = [];
cfg.tvec = lfp.tvec;
cfg.SWR = spwr_iv;
cfg.MethodSelection = [1 0 0 0]; %use: help getCandSeq
CAND_iv = getCandSeq(cfg,S_pc,pos);

%% Plot candidate events
cfg = [];
cfg.evt = CAND_iv;
cfg.windowSize = 1;
MultiRaster(cfg,S_pc)

%% Score candidate events
cfg = [];
cfg.display = 0;
scores = scoreCandSeq(cfg,CAND_iv,S_pc);

%% Plot significant events

% add scores to usr field of candidate ivs. N.B.- this can be done in scoreCandSeq() but
% for now we have a "scores" output which we can append back onto the candidate ivs.
datas = struct2cell(scores);
label = fieldnames(scores);
for f = 1:length(label)
    CAND_iv.usr(f).data = datas{f}'; %transpose to make it a column vector
    CAND_iv.usr(f).label = label{f};
end

% The navigate function now includes usr field text in the title for each event so
% MultiRaster will show the replay scores for each event. 
cfg = [];
cfg.lfp = lfp;
cfg.evt = CAND_iv;
MultiRaster(cfg,S_pc)



