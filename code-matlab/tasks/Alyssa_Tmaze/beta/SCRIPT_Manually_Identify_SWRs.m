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
% will open showing you the full length of data. The local field potential
% is plotted at the bottom of the window in black. Intervals where the ripple
% envelope exceeds the mean are plotted in red on the LFP, and this can
% serve as a hint to narrow down what could be a ripple.
% Press the 'b' key to go to the beginning of the data. Now, hit the 
% "Teleport" button on the far right of the window with the drop down menus
% set to '1' and 'beginning'. Now you should see a vertical bar in the
% middle of the window with a purple patch object to the left and a white
% background to the right. (Do not identify SWRs that are in the purple
% regions unless the SWR straddles the boundary.)
% To scroll through the data, move the figure window right and
% left using the arrow keys. When you have reached the end of the first
% segment of data you will see another vertical bar. You do not have to
% navigate to the next segment using the keyboard (that might take a while
% or some luck). Instead, select '2' and 'beginning' in the dropdown menus
% and hit "Teleport". This takes you to the beginning of the second
% segment. You can navigate between red LFP intervals using the 'a' key to
% go left and the 'd' key to go right.
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

% Set your current directory to the session you want to work with

% Select t files or tt files
cfg.whichSFiles = 'tt'; % 't' for MClust t files (isolated units), or 'tt' for MakeTTFiles tt files (all spikes from tetrode)
% 'none' for no spike files available

% Do you want the LFP filtered?
cfg.FilterLFP = 0; % If 1, filters LFP between 40-450 Hz before plotting; if 0, doesn't

cfg.Restrict = 0;

% If this section describes what you want, hit "Run" or press F5

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                      MAIN BODY OF SCRIPT                            %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load some things, check some things

[~,sessionID,~] = fileparts(pwd);

% Identify rat
ratID = SWRratID(sessionID);

% get LFP
LoadExpKeys
if isFrankRat(ratID)
    LFP = LoadFrankEEG([],ExpKeys.goodSWR{1});
    
elseif isvanderMeerRat(ratID) 
    please = [];
    please.fc = ExpKeys.goodSWR(1);
    LFP = LoadCSC(please);
    
else
    error('LFP loading procedure has not yet been implemented for this rat')
end

% Load Spikes
switch cfg.whichSFiles
    case 'tt'
        % Load all tetrode spikes
        if isvanderMeerRat(ratID)
            please = [];
            please.fc = FindFiles('*.tt');
            S = LoadSpikes(please);
        else
            error('Loading procedure for all recorded spikes has not yet been implemented for this rat')
        end
    case 't'
        % Load isolated units
        if isvanderMeerRat(ratID)
            please = []; % load
            please.load_questionable_cells = 1;
            S = LoadSpikes(please);
        elseif isFrankRat(ratID)
            please = [];
            S = LoadFrankSpikes(please,ExpKeys.epoch);
        else
            error('Loading procedure for isolated spikes has not yet been implemented for this rat')
       end
    case 'none'
        S = ts;
        S.t{1}(1,1) = LFP.tvec(1); % make a fake S because MultiRaster requires this as an input in order to work
        S.t{1}(2,1) = LFP.tvec(end);
    otherwise
        error('Unrecognized cfg.whichSFiles. Better check that spelling.')
end

% get restriction times
cfg.nSeconds = 720; % 12 minutes of data
tend = LFP.tvec(end);
tstart = tend - cfg.nSeconds;

if tstart < LFP.tvec(1)
    tstart = LFP.tvec(1);
end

%% Filter LFP

if cfg.FilterLFP
    please = [];
    please.f = [40 450]; % in Hz
    please.type = 'fdesign'; % doit4me
    LFP = FilterLFP(please,LFP);
end

%% Restrict S and LFP

if cfg.Restrict
    S = restrict(S,tstart,tend);
    LFP = restrict(LFP,tstart,tend);
end

%% Generate header to store info

hdr.host = getenv('COMPUTERNAME'); % the computer's name that ran the script
hdr.MATLAB_version = version; % because just in case
hdr.date = datestr(now);
hdr.sessionID = sessionID; % we actually care about this one
hdr.LFP = LFP.label; % we actually care about this one
hdr.segments = iv(tstart,tend);
hdr.cfg = cfg;

%% CALL MANUALIV

please = [];
  
% Decide whether user is resuming from a previous session
fn = FindFiles(['*manualIV.mat']);
if length(fn) > 1
    error('More than one file matching *%s.mat was found.',cfg.fn)
elseif isempty(fn)
    % do nothing
else
    IVann = loadpop(fn{1});
    please.resume = IVann;
    assert(isequal(IVann.hdr.segments,hdr.segments))
end
    
please.hdr = hdr;
please.mode = 'unfixed';
please.segments = hdr.segments;
ManualIV(please,S,LFP)
