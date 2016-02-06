%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%              PAPER: Generate Tmaze Candidate Events               %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the pool of candidates that will be used for further
% analysis; relies on the function GenCandidateEvents()

% SWR & MUA detection
% Speed & Theta thresholding
% NActive cells thresholding
% ResizeIV (interval expansion/contraction) optional

% aacarey May 2015, Dec 2015 edit for paper (from MASTER_Generate..)

clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% Which rats to use?
cfg.rats = TmazeRats;

% Do you want to keep a record of command window history in the files subfolder for the session?
cfg.writeDiary = 1; % 1 - yes, 0 - no

% Prepend to the filename:
cfg.prefix = ''; % ex: 'test_'

% Append to filename
cfg.suffix = 'PAPER'; % [date,''];

% Do you want to see what the internal functions have to say?
cfg.verbose = 1; %1 - yes, 0 - no (if 0, GenCandidateEvents still talks, but no one else)

% Do you want to load the 5's?
cfg.load_questionable_cells = 1;

% Which SWR detector do you want to use? 
cfg.SWRmethod = 'HT'; % 'AM' for amSWR (frequency content similarity), 'HT' for OldWizard (hilbert transform), 'TR' for photonic (transient detection), or 'none' to skip

% Which MUA detector do you want to use?
cfg.MUAmethod = 'none'; % 'AM' for amMUA, or 'none' for skip MUA detection

cfg.weightby = 'amplitude'; % this applies to 'AM' amSWR and 'TR' photonic, but not 'HT' OldWizard

% If using amSWR, how fast do you want it to go?
cfg.stepSize = 4; % this applies to 'AM' SWR method only

% Threshold for making intervals
cfg.DetectorThreshold = 3; % the threshold you want for generating IV data

% How is the threshold done? 'raw','zscore'
cfg.ThreshMethod = 'zscore';

% Minimum duration of intervals
cfg.mindur = 0.02; % in seconds

% Disclude events when the rat was moving faster than this
cfg.SpeedLimit = 10; % ( suggested 10 ) pixels per second, if [] no speed thresholding

% Disclude events with high theta power
cfg.ThetaThreshold = 2; % ( suggested 2 ) power std above mean, if [] no theta limiting

% Minimum number of active cells for the event to be kept
cfg.minCells = 5;

% Amount to add to the interval (catch borderline missed spikes)
cfg.expandIV = [0 0]; % in seconds (see ResizeIV)

cfg.allowOverlap = 0; % don't allow the expanded intervals to overlap one another

%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.gen);
    % and requisites are already checked by overseer script, so don't do it
    % here.
else
    % verify requisites: detection is a long process, don't want it erroring partway through
    
    cfg_temp = [];
    cfg_temp.requireExpKeys = 1;
    cfg_temp.ExpKeysFields = {'goodSWR','goodTheta'};
    cfg_temp.requireMetadata = 1;
    cfg_temp.MetadataFields = {'SWRfreqs'};
    cfg_temp.requireVT = 1;
    cfg_temp.requireHSdetach = 1;
    if cfg.writeDiary, cfg_temp.requireFiles = 1; end
    
    if ~checkTmazeReqs(cfg_temp); return; end
end
%%

originalFolder = pwd;

cfg_temp = [];
cfg_temp.rats = cfg.rats;
fd = sort(getTmazeDataPath(cfg_temp)); % get all session directories

for iFD = 1:length(fd)
    cd(fd{iFD});
    
    [~,session,~] = fileparts(fd{iFD});
    
    savename = [cfg.prefix,session,'-candidates',cfg.suffix];
    
    if cfg.writeDiary % save command window text
        cd([fd{iFD},'\files'])
        diary([savename,'_text','.txt'])
        cd(fd{iFD})
    end
    
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~              Beginning new candidates run                   ~~~')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    fprintf('\nDATE of script run: %s\n',datestr(now));
    disp(' '); disp(['COMPUTERNAME: ',getenv('COMPUTERNAME')])
    
    % tell me what you're working on
    disp(['*** Working on ',session,' ***']);
    disp(' ');
    
    % generate candidates
    please.verbose = cfg.verbose;
    please.load_questionable_cells = cfg.load_questionable_cells;
    please.SWRmethod = cfg.SWRmethod;
    please.MUAmethod = cfg.MUAmethod;
    please.weightby = cfg.weightby;
    please.stepSize = cfg.stepSize;
    please.DetectorThreshold = cfg.DetectorThreshold;
    please.mindur = cfg.mindur;
    please.SpeedLimit = cfg.SpeedLimit;
    please.ThetaThreshold = cfg.ThetaThreshold;
    please.minCells = cfg.minCells;
    please.expandIV = cfg.expandIV;
    please.allowOverlap = cfg.allowOverlap;
    
    evt = GenCandidateEvents(please);
    
    if cfg.writeDiary, diary off, end
    
    % save candidates
    save([savename,'.mat'],'evt');
end

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                  End of candidates run                      ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

cd(originalFolder)