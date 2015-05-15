function evt = GenCandidateEvents(cfg_in)
%GENCANDIDATEEVENTS generate pool of candidate replay events
% 
% Includes: SWR detection, MUA detection, speed thresholding, theta
% thresholding, number of active cells thresholding, and interval expansion. 
%
% REQUIRED FIELDS
%     ExpKeys.goodSWR
%     ExpKeys.goodTheta
%     metadata.SWRfreqs
%
% CONFIG OPTIONS
%     open the function and see the section called "parse cfg parameters"
%
% ACarey Mar 2015, for T-maze candidate events

%% Tell me what you're doing

tic

cprintf(-[0 0 1],'GenCandidateEvents: Looking for candidate replay events:');
disp(' ');
%% parse cfg parameters

cfg_def.verbose = 0; % i don't want to see what the internal functions have to say 
cfg_def.load_questionable_cells = 1;
cfg_def.weightby = 'amplitude';
cfg_def.DetectorThreshold = 4; % the threshold you want precand to use
cfg_def.mindur = 0.02; % in seconds, the minumum duration for detected events to be kept
cfg_def.SpeedLimit = 10; % pixels per second
cfg_def.ThetaThreshold = 2; % power, std above mean
cfg_def.minCells = 5; % minimum number of active cells for the event to be kept
cfg_def.expandIV = [0 0]; % amount to add to the interval (catch borderline missed spikes)
cfg_def.allowOverlap = 0; % don't allow the expanded intervals to overlap one another

cfg = ProcessConfig2(cfg_def,cfg_in);

%% Load some data
disp(' ')
disp('Loading data')

LoadExpKeys
LoadMetadata % for freqs

cfg_temp.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(cfg_temp);

cfg_temp = []; cfg_temp.useClustersFile = 0; cfg_temp.load_questionable_cells = cfg.load_questionable_cells;
S = LoadSpikes(cfg_temp);

cfg_temp = []; cfg_temp.fc = ExpKeys.goodTheta(1);
lfp_theta = LoadCSC(cfg_temp);

pos = LoadPos([]);

%% SWR score
disp(' ')
disp('Working on SWR detection......')

cfg_temp = []; cfg_temp.verbose = cfg.verbose; cfg_temp.weightby = cfg.weightby;
[SWR,~,~] = amSWR(cfg_temp,metadata.SWRfreqs,csc);

%% MUA score
disp(' ')
disp('Working on MUA detection.....')

cfg_temp = []; 
cfg_temp.verbose = cfg.verbose;
[MUA,~,~] = amMUA(cfg_temp,S,csc.tvec);

%% Threshold the detector
disp(' ')
disp('Compiling and thresholding detection scores....')

cfg_temp =[];
cfg_temp.threshold = cfg.DetectorThreshold;
cfg_temp.verbose = cfg.verbose; cfg_temp.mindur = cfg.mindur;
evt = precand(cfg_temp,csc.tvec,SWR,MUA,S);

%% Speed thresholding
disp(' ')
disp('Limiting output based on speed...')

cfg_temp = [];
spd = getLinSpd(cfg_temp,pos);

cfg_temp.threshold = cfg.SpeedLimit; cfg_temp.dcn = '<'; cfg_temp.method = 'raw';
%cfg_temp.minlen = 0; 
low_spd_iv = TSDtoIV(cfg_temp,spd);

evt = restrict(evt,low_spd_iv);
message = [num2str(length(evt.tstart)),' speed-limited candidates found.'];
disp(message)

%% Theta power thresholding
% in case there's theta even when the rat is moving very slowly or is
% stationary?
disp(' ')
disp('Limiting output based on theta power..')

% remove events during high theta
% 4th order, passband [6 10]
cfg_temp = []; cfg_temp.f = [6 10]; cfg_temp.order = 4; cfg_temp.display_filter = 0; cfg_temp.type = 'fdesign';
lfp_theta_filtered = FilterLFP(cfg_temp,lfp_theta);

% create tsd with theta power (lfp_theta is the output from LoadCSC)
tpow = LFPpower([],lfp_theta_filtered);

% detect intervals with "low theta" (z-score below 2)
cfg_temp = []; cfg_temp.method = 'zscore'; cfg_temp.threshold = cfg.ThetaThreshold; cfg_temp.dcn =  '<';
low_theta_iv = TSDtoIV(cfg_temp,tpow);

% restrict candidate events to only those inside low theta intervals
evt_test = restrict(evt,low_theta_iv);

disp([num2str(length(evt_test.tstart)),' theta-limited candidates found.'])

%% Number of active cells thresholding (this is done last because it's slower than speed and theta thresholding)
disp(' ')
disp('Limiting output based on number of active cells.')

iv_in = iv(evt.tstart,evt.tend);
activeCellsIV = AddNActiveCellsIV([],iv_in,S);
evt.nActiveCells = activeCellsIV.usr.data;
exclude = arrayfun(@(x) x < cfg.minCells,evt.nActiveCells);
evt.tstart(exclude) = [];
evt.tend(exclude) = [];
evt.nActiveCells(exclude) = [];

%% Expand intervals to catch missed spikes
disp(' ')
disp('Expanding intervals of remaining events')

cfg_temp = [];
cfg_temp.d = cfg.expandIV;
evt = addIV(cfg_temp,evt);

%% Keep a record of parameters

parameters.DetectorThreshold = cfg.DetectorThreshold;
parameters.mindur = cfg.mindur; % in seconds, the minumum duration for detected events to be kept
parameters.SpeedLimit = cfg.SpeedLimit; % pixels per second
parameters.ThetaThreshold = cfg.ThetaThreshold; % power, std above mean
parameters.minCells = cfg.minCells; % minimum number of active cells for the event to be kept
parameters.expandIV = cfg.expandIV; % amount to add to the interval (catch borderline missed spikes)
parameters.allowOverlap = cfg.allowOverlap; 
evt.parameters = struct('parameters',parameters,'SWRfreak',SWR.parameters.SWRfreak,'amSWR',SWR.parameters,'amMUA',MUA.parameters);

% tell me how many you found
disp(['Finished detection: ',num2str(length(evt.tstart)),' candidates found.'])

toc
end

