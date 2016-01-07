%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%               Annotate Manually Identified Intervals                %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ~~~~~~~~~   INFORMATION  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Summary: You will rate or delete the intervals you previously identified.
%
% AnnotateIV
% This script relies on the function AnnotateIV, which has a basic ui that
% allows the user to move through intervals previously identified in neural
% data and add annotations. It is very similar to ducktrap.
%
% Rating system (colour descriptions are for 'root' color axis scaling)
%
% 1  -  DEEP RED peak between 140-250 Hz
%    -  "If this isn't a SWR, then I surrender my rights as a human."
%    -  Every expert in the field would identify this as a SWR.
%    -  Ripple is well formed with many high amplitude ripples and has a
%       profile shaped like a "bundle" or "spindle" (i.e. with peak ripple
%       height near the middle of the event), and you definitely don't need
%       a spectrogram to tell it's a SWR.
%
% 2  -  RED/ORANGE (medium red or orange) peak between 140-250 Hz
%    -  "Even amateurs would recognize this as a SWR after little practice."
%    -  Ripple is well formed with many high amplitude ripples and has a
%       profile shaped like a "bundle" or "spindle" (i.e. with peak ripple
%       height near the middle of the event), and you definitely don't need
%       a spectrogram to tell it's a SWR.
%
% 3  -  YELLOW or GREEN peak between 140-250 Hz
%    -  "Most people would identify this as a moderate SWR."
%    -  Ripple is distinguishable from surrounding LFP but may be of short 
%       duration or moderate amplitude. It is recognizable as a SWR
%       outside of the context of surrounding LFP and MUA, and you probably
%       don't need a spectrogram to identify it as a SWR.
%
% 4  -  CYAN peak between 140-250 Hz
%    -  "Maybe or mostly a SWR, or possibly just a high-frequency event (HFE)."
%    -  Ripple amplitude is slightly higher than surrounding LFP.
%    -  Event may be of short duration or long and somewhat poorly formed.
%       Still recognizable as a SWR to a trained eye outside of the context 
%       of MUA, though it may require a spectrogram to distinguish it from
%       a 5.
%
% 5  -  LIGHT BLUE peak between 140-250 Hz
%    -  "This is a high frequency event (HFE) but probably doesn't make the
%       cutoff as a true SWR according to most people. Deleting it would
%       give me some amount of anxiety."
%    -  Ripple amplitude is similar or slightly higher than surrounding LFP.
%    -  Event may be of short duration or long and poorly formed. It may
%       not be recognizable as a SWR/HFE without the context of surrounding
%       LFP and MUA or the aid of a spectrogram.
%    -  Despite the above uncertainty as to whether this is a SWR, you ARE
%       fairly certain that this does not represent the spiking activity of
%       a single cell on this one tetrode. (i.e. given a similar quality
%       tetrode from another dCA1 recording site, you would still expect to
%       find frequency content similar to what you have observed here).
%
% delete  -  MEDIUM/DARK BLUE peak between 140-250 Hz in 'root' mode, or peak
%            disappears in 'raw' mode.
%         -  "What was I thinking?" or "I really don't know, it barely
%            makes the cutoff but keeping it gives me anxiety."
%         -  There are so many events like this that you think it's 
%            unreasonable to keep it, or
%         -  There is a peak in the proper range, but you can see visible
%            spikes in the LFP that distinctly correspond to spikes on the
%            corresponding tetrode in the rasterplot. There really isn't
%            much MUA compared to events labeled 1-5, or
%         -  No peak in the SWR range. Maybe you made a mistake before, or
%            seeing a lot of SWRs has helped you better recognize what a
%            SWR really looks like. 
%
% What if my interval matches one of the above descriptions perfectly
% except that the peak frequency range is between 100-140 Hz?
%    -  If there is another peak within the proper range that is not the
%       result of a single cell's spikes, assign your ratings based on the
%       second peak
%    -  If the only peak is between 100-140 Hz, dock a rating point or call
%       it a 4 or 5 (maybe a SWR)
%
% What if my interval is technically a 4 according to the rating system,
% but I am very confident that it's truly a SWR?
%    - Give it a rating of 3
%
% What if my interval has one really nice peak between 50-100 HZ?
%    -  You've found a really nice gamma event.
%
%
% ~~~~~~~~~   INSTRUCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% See config options in the section below titled "WHAT DO YOU WANT THIS SCRIPT
% TO DO?"
% Run the script.
% A figure window will open with various UI control buttons. One of your
% previously identified intervals will be plotted at the center of the
% window between vertical red lines. This is the interval you will provide
% a rating for. (You may see neighbouring intervals plotted in red in the
% LFP.) 
% On the right side of the window are the buttons "Previous" and "Next"
% between which are an edit box for you to type your rating and a checkbox
% to select whether you want your interval to be deleted. If you want to
% modify the interval boundaries, click the "Redraw" button. You will be
% prompted to click on the beginning and then the end of the interval. If
% you want to exit the redraw process before completion hit "Cancel".
% Each time you rate or select to delete an interval, the Progress updates
% at the top right of the window.
% If you want to display only the only the intervals you have not annotated
% (or selected for deletion), choose "Show unseen" from the drop down menu.
% To see all the intervals choose "Show all".
% 
% You can use the navigate function to scroll through the figure window,
% though this is discouraged: AnnotateIV wants you to annotate the interval
% it has shown you, and if you navigate away it is still considering that
% same interval but is giving you an opportunity to look around first. Make
% sure you go back to that interval by hitting [Previous then Next] or
% [Next then Previous] to see the current interval again. In general,
% navigation should not be necessary, so no "find current event" button was
% created.
% Note that clicking some UI controls removes focus from the axes and 
% navigate will not respond unless you return focus to the axes by clicking
% on the plot.
%
% In addition to the keyboard commands used by navigate, there are some
% that are specific to AnnotateIV:
%
%  s - spectrogram
%  n - go to the 'next' event (not the same as navigate's 'd')
%  p - go to the 'previous' event (not the same as navigate's 'a')
%  t - 'type', go to the annotate box (when in the annotation box, hit Enter to
%      return focus to the figure window and continue using keyboard commands)
%      If you're wondering whether AnnotateIV kept what you typed, the
%      answer is probably yes: as long as there was text entered in the box
%      at any time, this text was kept whether Enter was pressed or not.
%      You can check by going back.
%  r - 'remove', delete this event. You can change your mind about deleting
%      an event at any time as long as the figure window remains open.
%      Events are deleted upon saving: If you close the figure, the intervals 
%      are exactly as they were the last time you saved them.
%
%  There are no keyboard shortcuts for Redraw or switching to unseen
%  intervals.
% 
% Save your progress frequently.
%
% aacarey Dec 2015

clear

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% Set your current directory to the session you want to work with.

% What is the filename for your manually identified intervals?
cfg.fn = 'manualIV'; % unique string identifier for the file you containing iv data you want to annotate

% Do you want to be shown the intervals in random or chronological order?
cfg.mode = 'random'; % 'chrono' or 'random'

% Do you want to filter the LFP (this helps you see the shape of the SWR
% by removing underlying oscillations or w/e)
% If empty [], uses the settings saved from manual identification
% If 1, filters LFP
% If 0, does not filter the LFP
cfg.FilterLFP = 1; % [], 0, or 1

% If this section describes what you want, hit "Run" or press F5

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                      MAIN BODY OF SCRIPT                            %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load some things, check some things

fn = FindFiles(['*',cfg.fn,'.mat']);
if length(fn) > 1
    error('More than one file matching *%s.mat was found.',cfg.fn)
elseif isempty(fn)
    error('Zero files matching *%s.mat were found.',cfg.fn)
else
    evt = loadpop(fn{1});
end

% Load Spikes
switch evt.hdr.cfg.whichSFiles
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
        error('Unrecognized cfg.whichSFiles. This should not even happen.')
end

please = [];
please.fc = evt.hdr.CSC;
CSC = LoadCSC(please);

%% Filter LFP

please = [];

switch cfg.FilterLFP
    case []
        please.FilterLFP = evt.hdr.cfg.FilterLFP;
    case 0
        please.FilterLFP = 0;
    case 1
        please.FilterLFP = 1;
    otherwise
        error('Unrecognized cfg.FilterLFP option specified')
end

if please.FilterLFP
    please = [];
    please.f = [40 450]; % in Hz
    please.type = 'fdesign'; % doit4me
    CSC = FilterLFP(please,CSC);
end

%% run AnnotateIV

please = [];
please.mode = cfg.mode;

AnnotateIV(please,evt,S,CSC)
