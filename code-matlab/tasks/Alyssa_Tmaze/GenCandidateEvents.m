function evt = GenCandidateEvents(cfg_in)
%GENCANDIDATEEVENTS generate pool of candidate replay events
% 
% Includes: SWR detection, MUA detection, speed thresholding, theta
% thresholding, number of active cells thresholding, and interval expansion. 
%
% REQUIRED FIELDS
%     ExpKeys.goodSWR
%     ExpKeys.goodTheta
%     metadata.SWRfreqs (if using amSWR)
%
% CONFIG OPTIONS
%   cfg.verbose = 0; % i don't want to see what the internal functions have to say 
%   cfg.load_questionable_cells = 1;
%   cfg.SWRmethod = 'AM'; % 'AM' for amSWR (frequency content similarity), 
%       'HT' for OldWizard (hilbert transform), 'TR' for photonic 
%       (transient detection), or 'none' to skip
%   cfg.MUAmethod = 'AM'; % 'AM' for amMUA, or 'none' for skip MUA detection
%   cfg.weightby = 'amplitude'; % 'power' or 'amplitude', this applies to 
%       'AM' amSWR and 'TR' photonic, but not 'HT' OldWizard
%   cfg.stepSize = 1; % if using amSWR, this controls the speed of
%       detection.
%   cfg.ThreshMethod = 'zscore'; % 'zscore' or 'raw' how you want the
%       detector's TSD to be thresholded
%   cfg.DetectorThreshold = 3; % the threshold you want to use for
%       thresholding the detector to generate intervals
%   cfg.mindur = 0.02; % in seconds, the minumum duration for detected 
%       events to be kept
%   cfg.SpeedLimit = 10; % max speed allowed in pixels/sec. If empty [],
%       this step is skipped.
%   cfg.ThetaThreshold = 2; % allowable max theta power, std above mean, If
%       empty [], this step is skipped
%   cfg.minCells = 5; % minimum number of active cells for the event to be 
%       kept. if this is empty [], this step is skipped
%   cfg.expandIV = [0 0]; % amount to add to the interval (catch borderline 
%       missed spikes). This is the cfg.d input to ResizeIV
%   cfg.allowOverlap = 0; % don't allow the expanded intervals to overlap one another
%
% ACarey Mar 2015, for T-maze candidate events
% aacarey edit Jan 2016

%% Tell me what you're doing

tic
mfun = mfilename;
[~,sessionID,~] = fileparts(pwd);
disp(['***',mfun,': Looking for candidate events:']);

%% parse cfg parameters

cfg_def.verbose = 0; % i don't want to see what the internal functions have to say 
cfg_def.load_questionable_cells = 1;
cfg_def.SWRmethod = 'AM'; % 'AM' for amSWR (frequency content similarity), 'HT' for OldWizard (hilbert transform), 'TR' for photonic (transient detection), or 'none' to skip
cfg_def.MUAmethod = 'AM'; % 'AM' for amMUA, or 'none' for skip MUA detection
cfg_def.weightby = 'amplitude'; % this applies to 'AM' amSWR and 'TR' photonic, but not 'HT' OldWizard
cfg_def.stepSize = 1;
cfg_def.ThreshMethod = 'zscore';
cfg_def.DetectorThreshold = 3; % the threshold you want for generating IV data
%if strcmp(sessionID,'R042-2013-08-17') || strcmp(sessionID,'R044-2013-12-23')
    %cfg_def.DetectorThreshold = 2.5;
%end
cfg_def.mindur = 0.02; % in seconds, the minumum duration for detected events to be kept
cfg_def.SpeedLimit = 10; % pixels per second
cfg_def.ThetaThreshold = 2; % power, std above mean
cfg_def.minCells = 5; % minimum number of active cells for the event to be kept. if this is empty [], this step is skipped
cfg_def.expandIV = [0 0]; % amount to add to the interval (catch borderline missed spikes)
cfg_def.allowOverlap = 0; % don't allow the expanded intervals to overlap one another

cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if strcmp(cfg.SWRmethod,'none') && strcmp(cfg.MUAmethod,'none')
    error('cfg.SWRmethod and cfg.MUAmethod are set to none. You have to select at least one.')
end

%% Load some data
disp(' ')
disp('***Loading data')

LoadExpKeys
LoadMetadata % for freqs

cfg_temp = []; cfg_temp.getRatings = 0; cfg_temp.load_questionable_cells = cfg.load_questionable_cells; cfg_temp.verbose = cfg.verbose;
S = LoadSpikes(cfg_temp);

if ~ isempty(cfg.SpeedLimit)
    cfg_temp = []; cfg_temp.verbose = cfg.verbose;
    pos = LoadPos(cfg_temp);
end

cfg_temp = []; cfg_temp.verbose = cfg.verbose;
cfg_temp.fc = ExpKeys.goodSWR(1);
CSC = LoadCSC(cfg_temp);

if ~isempty(cfg.ThetaThreshold) % load a theta CSC
    cfg_temp = []; cfg_temp.fc = ExpKeys.goodTheta(1); cfg_temp.verbose = cfg.verbose;
    lfp_theta = LoadCSC(cfg_temp);
end

%% SWR score
if ~strcmp(cfg.SWRmethod,'none'); disp(' '); disp('***Working on SWR detection......'); end

switch cfg.SWRmethod
    case 'AM'
        cfg_temp = []; cfg_temp.verbose = cfg.verbose; cfg_temp.weightby = cfg.weightby; cfg_temp.stepSize = cfg.stepSize;
        [SWR,~,~] = amSWR(cfg_temp,metadata.SWRfreqs,CSC);
        
    case 'HT'
        cfg.stepSize = [];
        cfg.weightby = [];
        cfg_temp = []; cfg_temp.verbose = cfg.verbose; cfg_temp.rippleband = [140 250]; cfg_temp.smooth = 1; cfg_temp.kernel = [];
        SWR = OldWizard(cfg_temp,CSC);
        
    case 'TR'
        cfg.stepSize = [];
        
        cfg_temp =[]; cfg_temp.type = 'fdesign'; cfg_temp.f = [140 250];
        CSCr = filterLFP(cfg_temp,CSC);
        
        cfg_temp.weightby = cfg.weightby; cfg_temp.kernel = 'gauss';
        SWR = photonic(cfg_temp,CSCr);
        
    case 'none'
        cfg.stepSize = [];
        cfg.weightby = [];
        SWR = [];
end

%% MUA score
if ~strcmp(cfg.MUAmethod,'none'); disp(' '); disp('***Working on MUA detection.....'); end

switch cfg.MUAmethod
    case 'AM'
        cfg_temp = [];
        cfg_temp.verbose = cfg.verbose;
        [MUA,~,~] = amMUA(cfg_temp,S,CSC.tvec);
        
    case 'none'
        MUA = [];
end
%% Combine the detectors
disp(' ')
disp('***Compiling and thresholding detection scores....')

if ~isempty(SWR)
    SWR = rescmean(SWR,1); % this was a step in precand()
elseif isempty(SWR) && ~isempty(MUA)
    SWR = MUA; % for getting geometric mean...if SWR == MUA , then gm(SWR,MUA) = MUA
end

if ~isempty(MUA)
    %MUA = rescmean(MUA,1); % rescaling of MUA was not performed in precand()
elseif isempty(MUA) && ~isempty(SWR)
    MUA = SWR;
end

cfg_temp =[];
cfg_temp.method = 'geometricmean'; cfg_temp.verbose = cfg.verbose;
score = MergeTSD(cfg_temp,SWR,MUA);

score = rescmean(score,0.5); % this was a step in precand(), doing same to preserve

% threshold
cfg_temp = [];
cfg_temp.threshold = cfg.DetectorThreshold; cfg_temp.verbose = cfg.verbose; cfg_temp.method = cfg.ThreshMethod;
cfg_temp.operation = '>';
evt = TSDtoIV2(cfg_temp,score);

% remove short intervals
cfg_temp = []; cfg_temp.verbose = cfg.verbose; cfg_temp.mindur = cfg.mindur;
evt = RemoveIV(cfg_temp,evt);

evt.data = score.data; evt.tvec = score.tvec;
% cfg_temp =[];
% cfg_temp.threshold = cfg.DetectorThreshold;
% cfg_temp.verbose = cfg.verbose; cfg_temp.mindur = cfg.mindur;
% evt = precand(cfg_temp,CSC.tvec,SWR,MUA,S);



%% Speed thresholding

if ~isempty(cfg.SpeedLimit)
    disp(' ')
    disp('***Limiting output based on speed...')
    
    cfg_temp = []; cfg_temp.verbose = cfg.verbose;
    spd = getLinSpd(cfg_temp,pos);
    
    cfg_temp.threshold = cfg.SpeedLimit; cfg_temp.dcn = '<'; cfg_temp.method = 'raw'; cfg_temp.verbose = cfg.verbose;
    %cfg_temp.minlen = 0;
    low_spd_iv = TSDtoIV(cfg_temp,spd);
    
    evt = restrict(evt,low_spd_iv);
    disp(['***',num2str(length(evt.tstart)),' speed-limited candidates found.'])
end

%% Theta power thresholding
% in case there's theta even when the rat is moving very slowly or is
% stationary?

if ~isempty(cfg.ThetaThreshold)
    disp(' ')
    disp('***Limiting output based on theta power..')
    
    % remove events during high theta
    % 4th order, passband [6 10]
    cfg_temp = []; cfg_temp.f = [6 10]; cfg_temp.order = 4; cfg_temp.display_filter = 0; cfg_temp.type = 'fdesign'; cfg_temp.verbose = cfg.verbose;
    lfp_theta_filtered = FilterLFP(cfg_temp,lfp_theta);
    
    % create tsd with theta power (lfp_theta is the output from LoadCSC)
    cfg_temp = []; cfg_temp.verbose = cfg.verbose;
    tpow = LFPpower(cfg_temp,lfp_theta_filtered);
    
    % detect intervals with "low theta" (z-score below 2)
    cfg_temp = []; cfg_temp.method = 'zscore'; cfg_temp.threshold = cfg.ThetaThreshold; cfg_temp.dcn =  '<'; cfg_temp.verbose = cfg.verbose;
    low_theta_iv = TSDtoIV(cfg_temp,tpow);
    
    % restrict candidate events to only those inside low theta intervals
    evt = restrict(evt,low_theta_iv);
    
    disp(['***',num2str(length(evt.tstart)),' theta-limited candidates found.'])
end

%% if there was HS detachment, we need to restrict evt

% 1. restrict based on HS_detach_times. note this is basically unnecessary
% because there are no spikes during these times (ntt file restricted),
% which means the MUA detector will zero out this region anyway
if strcmp(sessionID,'R044-2013-12-21') || strcmp(sessionID,'R044-2013-12-22')
    disp('***Excluding HS detach times')
    load(FindFile('*HS_detach_times.mat'))
    temp_iv = iv(t_start/10^6,t_end/10^6);
    %temp_iv.tstart(1) = lfp_theta.tvec(1);
    %temp_iv.tend(end) = lfp_theta.tvec(end);
    evt = restrict(evt,temp_iv);
end
% 2. restrict based on metadata.detachIV
if isfield(metadata,'detachIV')
    disp('***Restricting evt to metadata.detachIV')
    evt = restrict(evt,metadata.detachIV);
end

%% Number of active cells thresholding (this is done last because it's slower than speed and theta thresholding)
disp(' ')
cfg_temp =[]; cfg_temp.verbose = cfg.verbose;
evt = AddNActiveCellsIV(cfg_temp,evt,S);
if ~isempty(cfg.minCells)
    
    disp('***Limiting output based on number of active cells.')
    cfg_temp = []; cfg_temp.operation = '>='; cfg_temp.threshold = cfg.minCells;
    evt = SelectIV(cfg_temp,evt,'nActiveCells');
end

%% Expand intervals to catch missed spikes
if abs(prod(cfg.expandIV))> 0
    disp(' ')
    disp('***Expanding intervals of remaining events')
    
    cfg_temp = [];
    cfg_temp.d = cfg.expandIV; cfg_temp.allowOverlap = cfg.allowOverlap;
    evt = ResizeIV(cfg_temp,evt);
    if ~cfg.allowOverlap
        evt = MergeIV([],evt);
    end
end

% tell me how many you found
disp(' '); disp(['***Finished detection: ',num2str(length(evt.tstart)),' candidates found.']); disp(' ')

evt = History(evt,mfun,cfg);

toc
end

