%% Run event detection on R050-2014-04-02
% example

cd('D:\data\R050\R050-2014-04-02')

% get required inputs

fn = FindFiles('*metadata.mat');

load(fn{1})

cfg.useClustersFile = 0;
S = LoadSpikes(cfg);

cfg.fc = {'R050-2014-04-02-CSC07a.ncs'}; 
csc = LoadCSC(cfg);

%% SWR detector
% this does 2 ffts, so takes longer than original detectSWR_EG
% 10 mins

[SWR,~,~] = amSWR([],metadata.SWRfreqs,csc);
% replace tildas to get the two swr1 and swr2 that were merged to get SWR
% by default, these use 40 ms (precision) and 80 ms (reduce false
% positives) windows

%% MUA detector
% counts each spike as a discrete Dirac (as original detectMUA did), but
% has additional adaptive thresholding to decrease the frequency and breadth
% false positives.
% 20 seconds

[MUA,~,~] = amMUA([],S,csc.tvec);

% replace tildas to get "raw" and "noise" scores that were merged to get
% MUA
% raw is same as original detectMUA output. noise is the curve that does
% adaptive thresholding. 

%% combine and threshold
% 0.2 seconds if minCells not specified
% 14 seconds if specified

cfg.minCells = 5; % if want to exclude evt with too few active cells
cfg.threshold = 8;
evt = precand(cfg,csc.tvec,SWR,MUA,S); 

%% Final analysis function combines all three of the above functions; use when optimal cfg parameters are known:

 %[evt2,SWR2,MUA2] = detector(cfg,metadata.SWRfreqs,csc,S); 


%% Plot results

% set up MR
cfg_mr.evt = evt;
cfg_mr.lfp = csc;
cfg_mr.lfpHeight = 30;

% set up additional visualization aid (makes navigating slower when more
% plotted)

evt.name = 'interval'; % to show evt-inclusive spikes
evtscore = evt.data;

MultiRaster(cfg_mr,S);

sidekick(csc.tvec,evtscore,evt)

% view additional scores, useful for tweaking cfg parameters
%MultiRaster(cfg_mr,S);
%swrscore = SWR.data*3;
%muascore = MUA.data*3;
%sidekick(csc.tvec,evtscore,evt,swrscore,muascore)

