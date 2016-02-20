%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%        OVERSEER: Coordinate Tmaze Data Analysis Scripts             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Run this script to analyze data for the motivationl T-maze project.
%
%                    **************************
%                            IMPORTANT
%                    **************************
%
% If your promoted folders are not up to date you will likely run into issues.
%
% Your path must include the shared, Tmaze, and Replay_Analysis folders 
% from GitHub.
%
% Your current directory should be set up to have the following subfolders:
%    \code    (optional, but save a new version of this script here so it
%              stays with the data and visuals you made)
%    \data    (mandatory, for .mat and text files)
%    \visuals (mandatory, for figures)
%
%  Consider naming the main folder that holds these RunA and so on, or FINAL:
%   FINAL\code
%   FINAL\data
%   FINAL\visuals
%
%   Archive the folder that contains published results!
%
% Running this script with the current directory set to a fresh folder is
% recommended but not necessary if you keep track of what you have done
% already.
% 
%                    **************************
%
% This script coordinates the following Master scripts:
%   PAPER_Generate_Candidates
%   PAPER_Generate_inputData
%   PAPER_CoOccurrence
%   PAPER_CoOccurrence_ExportToR
%   PAPER_Plot_CoOccurrence
%   PAPER_Plot_CoOccurrence_Difference
%   PAPER_Collect_Behavior
%   PAPER_Collect_nEvents
%   PAPER_Collect_EventRate
%   PAPER_Collect_nUnitsRecorded
%   PAPER_Collect_nPlaceCells
%
% Each of these scripts can be run individually if the need arises, but be
% sure that the settings across scripts agree!
%
% Select master config options below. The master, CFG, is saved in the 
% \data folder as a record.
%
% aacarey Jan 2016

clear 

%% ~~~~~~~~~~         SELECT MASTER CONFIG OPTIONS           ~~~~~~~~~~~~~~

% GENERAL
% Do you want to use existing candidates?
CFG.SkipDetection = 0; % 1 - yes, 0 - no

% Do you want to use existing inputData?
CFG.SkipInputData = 0; % 1 - yes, 0 - no

% Do want to rerun CoOccurrence?
CFG.SkipCoOccurrence = 0; % 1 - yes, 0 - no

% CANDIDATES (see PAPER_Generate_Candidates for config options)
CFG.gen.suffix = 'PAPER'; % unique identifier for the candidates
CFG.gen.load_questionable_cells = 1;
CFG.gen.SWRmethod = 'HT';
CFG.gen.MUAmethod = 'none';
CFG.gen.stepSize = 4;
CFG.gen.DetectorThreshold = 3;
CFG.gen.ThreshMethod = 'zscore';
CFG.gen.SpeedLimit = 10;
CFG.gen.ThetaThreshold = 2;
CFG.gen.mindur = 0.02;
CFG.gen.minCells = []; % see also the same field in cooccurrence; you can do this here or at the cc step

% INPUTDATA (see PAPER_Generate_inputData for config options)
CFG.in.load_questionable_cells = 1;

% COOCCURRENCE (see PAPER_CoOccurrence for config options)
CFG.cc.whichCandidates = ['-candidates',CFG.gen.suffix];
% note that some scripts require all of {'prerecord','task','equalBehaviorITI','postrecord'}
CFG.cc.WhichEvents = {'prerecord','task','equalBehaviorITI','postrecord'};
CFG.cc.whichS = 'unique';
CFG.cc.min_cluster_quality = 5;
CFG.cc.minCells = 2;
CFG.cc.useMask = 1; % 1, final
CFG.cc.nShuffle = 10000; % use 10 000 for final run
CFG.cc.WhichP = {'p0','p4','p5'};

% PLOTTING
CFG.plot.colormode = 'inventory4'; % color scheme
CFG.plot.whichCandidates = CFG.cc.whichCandidates;

%~~~~~~~~~~~               END OF CONFIG OPTIONS               ~~~~~~~~~~~~

%% MAIN BODY OF SCRIPT

%~~~~~~~~ save some things as a record ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CFG.date = datestr(now);
CFG.host = getenv('COMPUTERNAME');
CFG.user = getenv('USERNAME');
CFG.MATLABversion = version;

%~~~~~~~~ make sure subfolders exist by trying to cd to them ~~~~~~~~~~~~~~

originalFolder = pwd;
cd([originalFolder,'\visuals'])
cd([originalFolder,'\data'])
cd(originalFolder)

CFG.mainFD = pwd;

if ~CFG.SkipDetection || ~CFG.SkipDetection
    
    %~~~~~~~~~~~ verify existence of some requisite data ~~~~~~~~~~~~~~~~~~
    cfg_temp = []; cfg_temp.checkall = 1;
    if ~checkTmazeReqs(cfg_temp); return; end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~ Run scripts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ~CFG.SkipDetection
    % SWR detection
    run('PAPER_Generate_Candidates')
end

if ~CFG.SkipInputData
    % Get tuning curves, place cells, etc...
    run('PAPER_Generate_inputData')
end

if ~CFG.SkipCoOccurrence
    % Perform co-occurrence analysis
    for iEvents = 1:length(CFG.cc.WhichEvents)
        % get set of co-occurrence results for this group of events
        
        CFG.cc.whichEvents = CFG.cc.WhichEvents{iEvents};
        
        run('PAPER_CoOccurrence')
        
        CFG.cc = rmfield(CFG.cc,'whichEvents');
    end
end

% Handle some co-occurrence data
for iP = 1:length(CFG.cc.WhichP)
    % change which p field we're working with (refer to CoOccurQ output)
    
    CFG.plot.whichP = CFG.cc.WhichP{iP};
    
    % Format data for R stats
    run('PAPER_CoOccurrence_ExportToR')
    
    % Plot raw co-occurrence results
    run('PAPER_Plot_CoOccurrence');
    
    % Plot processed co-occurrence results
    run('PAPER_Plot_CoOccurrence_Difference');
    
    CFG.plot = rmfield(CFG.plot,'whichP');
end

% Trial choice, and plot
run('PAPER_Collect_Behavior')

% Total SWR counts
run('PAPER_Collect_nEvents')

% SWRs per minute, plot
run('PAPER_Collect_EventRate')

% How many neurons were recorded from, LaTeX table
run('PAPER_Collect_nUnitsRecorded')

% How many place cells, LaTeX table
CFG.temp.threshold = CFG.cc.min_cluster_quality;
CFG.temp.operation = '<=';
run('PAPER_Collect_nPlaceCells')
CFG = rmfield(CFG,'temp');

cd([CFG.mainFD,'\data'])
save('MASTER_CONFIG','CFG') % in \data
cd(CFG.mainFD)

clearvars -except CFG

disp(datestr(now))
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~              END OF ANALYSIS RUN                 ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
