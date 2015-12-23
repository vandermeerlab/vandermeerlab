%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                   Identify SWRs Humanually                          %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This the slowest way of detecting SWRs, except perhaps for training a
% paramecium to do it for you instead. (Also note the very inclusive/PC 
% title of this script, unless you actually are a paramecium.)


% ~~~~~~~~~   INFORMATION  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Summary: You will click on the start and end times of SWRs. The
% available information includes LFP spikes, and a spectrogram option. The 
% interface is an interactive figure window.

% About ducktrap()
% This script relies on the function ducktrap (which may have been renamed
% to ManualIV). ducktrap has a basic UI that allows the user to scroll
% through neural data and identify regions of interest (such as SWRs). This
% script will use ducktrap's 'unfixed' mode, which means that you must
% click on the beginning and end of an interval (it's okay if you click on
% the end first: ducktrap will police your physics and fix the paradox).
% When you have identified an interval, it will appear in orange on the
% figure. To keep this interval hit "Enter" or "Return"; it will now appear
% in green. If you do not want to keep the interval just click anywhere else 
% on the figure and it will disappear. Don't get Enter-happy, because once
% you've hit Enter there is no way of undoing.
% Type "doc ducktrap" in the command window for more information.
%
% Tip: Watch where you click. If you click somewhere and navigate away, the
% plot is still there but off the viewing window. If you forget you've
% clicked once, then move and click again, there will be a plot appearing
% seemingly out of nowhere. Just click again to remove it.

% About navigate()
% Ducktrap calls on the keypress function 'navigate'. Navigate allows you
% to move the figure window using keyboard input. Navigate has detailed
% descriptions of the commands available to you. Especially useful ones are
% -> and <- which allow you to move the window right and left. 'w' allows you
% to change the viewing window size. If no interval data is plotted, the
% 'a' and 'd' keys move the window left and right in units of seconds.
% Adding shift and ctrl multiplies the distance moved by 10 and 50 seconds
% respectively. Type "doc navigate" in the command window for more
% information.

% ~~~~~~~~~   INSTRUCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% In this particular script, the data is divided into three segments: one for
% each of the prerecord, task, and postrecord. When you run this script for
% the first time, make sure that cfg.resumeSession = 0. A figure window
% will open showing you the full length of data. Hit 'f' on the keyboard to
% go to the first event (in this script, each event is one second of data).
% You should see a vertical bar in the center of the screen. To verify that
% this is the first segment of data, you can hit the "Teleport" button on 
% the far right of the window with the drop down menus set to '1' and 
% 'beginning'. To scroll through the data, move the figure window right and
% left using the arrow keys. When you have reached the end of the first
% segment of data you will see another vertical bar. You do not have to
% navigate to the next segment using the keyboard (that might take a while
% or some luck). Instead, select '2' and 'beginning' in the dropdown menus
% and hit "Teleport".
%
% Save your progress frequently because ducktrap is in beta and power
% outages don't care about you.
%
% If you want to continue from a previous session, change cfg.resumeSession
% = 1. If any config options aside from cfg.resumeSession and cfg.fn are
% different from the previous session, this script will error
% (intentionally). If your settings are the same, the figure will open to
% the last event you identified. Previously identified events are shown in
% blue.

% Please have your chocolate bars and epic music on hand.

% aacarey Dec 2015

clear

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% Set your current directory to the Tmaze session you want to work with

% T-maze specific recording segments to restrict data to
cfg.whichEpochs = {'prerecord','task','postrecord'}; %  {'prerecord','task','postrecord'} 

% How many seconds of data do you want to peruse? (this gets evenly divded between the chosen epochs)
cfg.nSeconds = 720; % 720 seconds = 12 minutes

% Select t files or tt files
cfg.whichSFiles = 'tt'; % 't' for MClust t files (isolated units), or 'tt' for MakeTTFiles tt files (all spikes from tetrode)

% Do you want the LFP filtered?
cfg.FilterLFP = 1; % If 1, filters LFP between 40-450 Hz bfore plotting; if 0, doesn't

% Do you want to restrict the data?
cfg.Restrict = 0; % If 1, restricts data to regions of interest, if 0 doesn't

% *** Would you like to resume from a previous session?
cfg.resumeSession = 0; % 1 if you want to resume, 0 if you're starting fresh
cfg.fn = 'manualIV'; % unique string identifier for the file you saved from a previous session

% If this section describes what you want, hit "Run" or press F5

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                      MAIN BODY OF SCRIPT                            %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load some things, check some things

% Load Spikes
switch cfg.whichSFiles
    case 'tt'
        % Load all tetrode spikes
        please = [];
        please.fc = FindFiles('*.tt');
        S = LoadSpikes(please);
    case 't'
       % Load isolated units
       please = []; % load 
       please.load_questionable_cells = 1;
       S = LoadSpikes(please);
    otherwise
        error('Unrecognized cfg.whichSFiles. Better check that spelling.')
end

% Just cuz
if any(~strcmp(cfg.whichEpochs,{'prerecord','task','postrecord'}))
    error('cfg.whichEpochs contains unrecognized recording segment')
end

% get CSC
LoadExpKeys
please = [];
please.fc = ExpKeys.goodSWR(1);
CSC = LoadCSC(please);

% get restriction times
nSecondsPerEpoch = cfg.nSeconds/(length(cfg.whichEpochs));

tstart = nan(size(cfg.whichEpochs));
tend = nan(size(cfg.whichEpochs));

for iEpoch = 1:length(cfg.whichEpochs)
    tstart(iEpoch) = ExpKeys.(cfg.whichEpochs{iEpoch})(2) - nSecondsPerEpoch;
    tend(iEpoch) = ExpKeys.(cfg.whichEpochs{iEpoch})(2);
end

%% Filter LFP

if cfg.FilterLFP
    please = [];
    please.f = [40 450]; % in Hz
    please.type = 'fdesign'; % doit4me
    CSC = FilterLFP(please,CSC);
end

%% Restrict S and LFP

if cfg.Restrict
    S = restrict(S,tstart,tend);
    CSC = restrict(CSC,tstart,tend);
end

%% Generate header to store info

hdr.host = getenv('COMPUTERNAME'); % the computer's name that ran the script
hdr.MATLAB_version = version; % because just in case
hdr.date = datestr(now);
[~,sessionID,~] = fileparts(pwd);
hdr.sessionID = sessionID; % we actually care about this one
hdr.CSC = CSC.label; % we actually care about this one
hdr.segments = iv(tstart,tend);
hdr.cfg = cfg;

%% CALL DUCKTRAP

please = [];

if cfg.resumeSession % load evt and check that settings agree
    
    fn = FindFiles(['*',cfg.fn,'.mat']);
    if length(fn) > 1
        error('More than one file matching *%s.mat was found.',cfg.fn)
    elseif isempty(fn)
        error('Zero files matching *%s.mat were found.',cfg.fn)
    else
        evt = loadpop(fn{1});
    end
    
%     hdr.verbose = 1;
%     ProcessConfig(evt.hdr.cfg,cfg); % hack to compare field names :P
%     hdr = rmfield(hdr,'verbose');
    
    % make sure session is the same, and tell user if not
    if ~isequal(evt.hdr.sessionID,sessionID)
        error('Reloaded events came from a different session: %s',sessionID)
    end
    
    % make sure CSC is the same, and tell user if not
    if ~isequal(evt.hdr.CSC,CSC.label)
        error('CSC from previous run and current run must be the same')
    end
    
    % bulk compare the whole hdr, if it's not the same tell the user, and
    % make them figure it out on their own
    temp_prevCFG = rmfield(evt.hdr.cfg,{'resumeSession','fn'});
    temp_currCFG = rmfield(hdr.cfg,{'resumeSession','fn'});
    
    if ~isequal(temp_prevCFG,temp_currCFG)
        error('headers for current run and previous run do not agree. see evt.hdr and compare cfg options')
    end
    
    please.resume = evt;
end

please.hdr = hdr;
please.mode = 'unfixed';
please.segments = hdr.segments;
ducktrap(please,S,CSC)
