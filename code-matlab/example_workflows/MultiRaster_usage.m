%% MULTRASTER EXAMPLE
% This script just shows how you can use the MultiRaster function. It takes a variety of
% event type objects (both iv and ts) as well as multiple tsd type data to create an
% interpertable plot. 


%% Load data and create event data
% This is taken mostly straight from the cosmo5 tutorial

clear all;
cd('D:\My_Documents\data\inProcess\R042-2013-08-18_recording'); %or wherever this folder is

% load a CSC (CA1 layer I think)
cfg = []; cfg.fc = {'CSC17.Ncs'};
lfp = LoadCSC(cfg);

% Load spikes
cfg = [];
cfg.useClustersFile = 0;
S = LoadSpikes(cfg);

% Load Position data
cfg = [];
pos_tsd = LoadPos(cfg);
cfg.realTrackDims = [167 185];
pos_test = LoadPos(cfg);

% Load events
cfg = [];
cfg.eventList = {'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0020).',...
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0080).',...
    'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0040).',...
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0004).',...
    'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0040).'};
cfg.eventLabel = {'center_pb','left_pb','right_pb','food_dispensed','water_dispensed'};
evt = LoadEvents(cfg);

% Load (or create) Coord data
load('D:\My_Documents\data\pfmodels\2014-11-04\coordL.mat');
load('D:\My_Documents\data\pfmodels\2014-11-04\coordR.mat');

% Resample Coord to get reasonable values
% NOTE - is this actually changing pixels to cm? it just looks like your binning the
% pixels into 180 bins so the value 
nBins = 180; % track is 180cm, aim for 1cm bins
CoordLrs(1,:) = interp1(1:size(CoordL,2),CoordL(1,:),linspace(1,size(CoordL,2),nBins),'linear');
CoordLrs(2,:) = interp1(1:size(CoordL,2),CoordL(2,:),linspace(1,size(CoordL,2),nBins),'linear');
CoordRrs(1,:) = interp1(1:size(CoordR,2),CoordR(1,:),linspace(1,size(CoordR,2),nBins),'linear');
CoordRrs(2,:) = interp1(1:size(CoordR,2),CoordR(2,:),linspace(1,size(CoordR,2),nBins),'linear');

cfg.Coord = CoordLrs;
posL = LinearizePos(cfg,pos_tsd);
cfg.Coord = CoordRrs;
posR = LinearizePos(cfg,pos_tsd);

clear CoordL CoordLrs CoordR CoordRrs nBins

% check left/right
load('times.mat');
leftruns = getd(evt,'water_dispensed');
for iRun = 1:length(run_start)
    if ~isempty(leftruns(leftruns > run_start(iRun) & leftruns < run_end(iRun)))
        evt.trial_ID(iRun) = 2; % water/right
    else
        evt.trial_ID(iRun) = 1; % food/left
    end 
end

evt.run_start = run_start;
evt.run_end = run_end;
evt.ped_start = ped_start*10^-6; % original data file did not convert to sec (s) so fix here
evt.ped_end = ped_end*10^-6;

left_trials = find(evt.trial_ID == 1);
evtL = evt;
evtL.trial_ID = evtL.trial_ID(left_trials);
evtL.run_start = evtL.run_start(left_trials);
evtL.run_end = evtL.run_end(left_trials);
evtL = restrict(evtL,evtL.run_start,evtL.run_end);
evtL = removeEmptyCells(evtL);

right_trials = find(evt.trial_ID == 2);
evtR = evt;
evtR.trial_ID = evtR.trial_ID(right_trials);
evtR.run_start = evtR.run_start(right_trials);
evtR.run_end = evtR.run_end(right_trials);
evtR = restrict(evtR,evtR.run_start,evtR.run_end);
evtR = removeEmptyCells(evtR);
 
Sl = restrict2(S,evtL.run_start,evtL.run_end);
Sr = restrict2(S,evtR.run_start,evtR.run_end);
posLrs = restrict(posL,evtL.run_start,evtL.run_end);
posRrs = restrict(posR,evtR.run_start,evtR.run_end);

clear left_trials leftruns right_trials prerecord_end prerecord_start t_end t_start 
clear ped_end ped_start run_end run_start iRun posL posR

%% TEST MULTIRASTER

% Create an iv for left runs
temp_iv = iv;
temp_iv.tstart = evtL.run_start;
temp_iv.tend = evtL.run_end;

% Create and split the lfp signal
start_time = lfp.tvec(1);
end_time = lfp.tvec(end);
splits = linspace(start_time,end_time,4);
lfpA = restrict2(lfp,splits(1),splits(2));
lfpB = restrict2(lfp,splits(1),splits(4));
lfpC = restrict2(lfp,splits(3),splits(4));

% Choose a plotting mode
% I've setup some example situations that one might come across when examining data. Try
% each and see! (Currently working on case 4 & 8).

plotmode = 7;
cfg = [];
cfg.windowSize = 10;
cfg.evtTimes = (evtL.run_start+evtL.run_end)./2;

tic;
switch plotmode
    case 1%original (just spikes)
        MultiRaster(cfg,S)
        
    case 2%ts events
        cfg.evt = evtL;
        cfg.legend = 'on';
        MultiRaster(cfg,S)
        
    case 3%iv events
        cfg.evt = temp_iv;
        MultiRaster(cfg,S)
        
    case 4%ts + iv events -- NOT WORKING
        error('Not working!')    
        
    case 5%lfp UPDATE -- can plot multiple lfps
        cfg.lfp(1) = lfpA;
        cfg.lfp(2) = lfpB;
        cfg.lfp(3) = lfpC;
        MultiRaster(cfg,S)
        
    case 6%lfp + ts
        cfg.lfp = lfp;
        cfg.evt = evtL;
        MultiRaster(cfg,S)
        
    case 7%lfp + iv
        cfg.lfp = lfp;
        cfg.evt = temp_iv;
        MultiRaster(cfg,S)
        
    case 8%lfp + ts + iv -- NOT WORKING
        error('Not working!')
end
toc;