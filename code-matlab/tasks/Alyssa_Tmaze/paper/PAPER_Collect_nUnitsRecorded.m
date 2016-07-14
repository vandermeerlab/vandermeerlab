%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%               COUNT TOTAL NUMBER OF UNITS RECORDED                  %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For Tmaze data
%
% This script loads .t and ._t files and counts all the units recorded per
% session for each rat. These are all the units recorded, except for those
% with very poor isolation or too few spikes to have been kept during the
% spike sorting process.
% Text is printed to the command window, including LaTeX code for a table,
% and saved to text file.
%
% Script assumes that each rat has exactly 6 session folders.
%
% aacarey Sept 2015 (poster version, two separate scripts)
%  --edit Jan 2016 (paper version, single script)
%  --mvdm ugly edit Jul 2016 (also count questionable units)

clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

cfg = [];

% Where to send the data and diary?
cfg.output_fd = [pwd,'\data']; 

% Which rats to include?
cfg.rats = TmazeRats;

% Should I load the 5's?
cfg.load_questionable_cells = 1;

% Would you like to save the text file?
cfg.writeFiles = 1; 

cfg.tablemode = 'default';

%% Info, data loading

originalFolder = pwd; % remember where we started

if cfg.writeFiles % save command window text
    cd(cfg.output_fd)
    diary(['nUnitsRecorded','_text','.txt'])
    cd(cfg.output_fd)
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~              COUNTING TOTAL NEURAL UNITS                    ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

fprintf('\nDATE of script run: %s\n',datestr(now));

%% Counting

% set up LoadSpikes' config
please = [];
please.load_questionable_cells = cfg.load_questionable_cells;

% Make something to collect each rat's session counts into. s1 is session 1
% and it will hold the number of units from that session.
empty = struct('s1',[],'s2',[],'s3',[],'s4',[],'s5',[],'s6',[],'total',[]);

TotalAllRats = 0; % grand number of units for experiment
TotalAllRatsQ = 0; % grand number of questionable units for experiment

% go through each rat and collect the data
for iRat = 1:length(cfg.rats);
    
    collect = empty; % make a copy of empty to put stuff into
    total = 0; % this will collect the number of units in total for a single rat
    totalQ = 0;
    
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['  Data for ',cfg.rats{iRat},':']); 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    % give me a list of directories where I can find session data:
    cfg_temp.rats = cfg.rats(iRat);
    fd = getTmazeDataPath(cfg_temp);  
    
    % go through each session and do the thing
    for iFD = 1:length(fd)
        
        cd(fd{iFD})
        [~,sessionID,~] = fileparts(fd{iFD});
        disp(' ')
        disp(['SessionID: ',sessionID,' (session #',num2str(iFD),')'])
        S = LoadSpikes(please);
        n = length(S.t);
        
        Snq = LoadSpikes([]); % no questionable cells
        nQ = n-length(Snq.t);
        
        % collect total cells recorded for all rats
        TotalAllRats = TotalAllRats + n;
        TotalAllRatsQ = TotalAllRatsQ + nQ;
        
        % collect total for each rat     
        total = total + n;
        totalQ = totalQ + nQ;
        
        % collect cells from each session
        collect.(['s',num2str(iFD)]).n = n; %  wow eh. 
        collect.(['s',num2str(iFD)]).nQ = nQ;
        % ^ tells matlab to assign the value to (for ex), empty.s1 for the first file directory
        % ['s',num2str(iFD)] if run, generates s1 or s2 and so on
        collect.(['s',num2str(iFD)]).sessionID = sessionID; % if needed for future data handling
        
    end
    
    collect.total = total; % the total cells recorded for the current rat
    collect.totalQ = totalQ; % the total cells recorded for the current rat
    nUnits.(cfg.rats{iRat}) = collect; % makes, for example, nUnits.R042 = collect;
    disp(' ')
    disp(['***Total units for ',cfg.rats{iRat},': ',num2str(total),' (',num2str(totalQ),')']); disp(' ')
end

nUnits.TotalAllRats = TotalAllRats;
nUnits.TotalAllRatsQ = TotalAllRatsQ;

disp(['Total units from all rats: ',num2str(nUnits.TotalAllRats),' (',num2str(nUnits.TotalAllRatsQ),')']); disp(' '); disp(' '); disp(' ')

%% print LaTeX code
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('      LATEX CODE FOLLOWS:'); 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'); disp(' ')


switch cfg.tablemode
    case 'default'
        disp('\begin{table}[h]')
        disp('\centering')
        disp('\caption[Total neural units]{\textmd{Total neural units for three rats across each of their six sessions.}}')
        disp('\small')
        disp('  \begin{tabular}{ l c c c c c c c c}')
        disp('    \hline')
        disp('    \textbf{Rat ID} & \textbf{Session 1}    &   \textbf{Session 2}    & \textbf{Session 3} & \textbf{Session 4}  & \textbf{Session 5}  & \textbf{Session 6} & \textbf{Total} \\')
		disp('    \hline')
		disp(['    \textbf{R042}    &     ',num2str(nUnits.R042.s1.n),' (',num2str(nUnits.R042.s1.nQ),')   &   ',num2str(nUnits.R042.s2.n),' (',num2str(nUnits.R042.s2.nQ),')  &  ',num2str(nUnits.R042.s3.n),' (',num2str(nUnits.R042.s3.nQ),')  &   ',num2str(nUnits.R042.s4.n),' (',num2str(nUnits.R042.s4.nQ),')  &   ',num2str(nUnits.R042.s5.n),' (',num2str(nUnits.R042.s5.nQ),')  &   ',num2str(nUnits.R042.s6.n),' (',num2str(nUnits.R042.s6.nQ),')  &\textbf{',num2str(nUnits.R042.total),' (',num2str(nUnits.R042.totalQ),')}\\ '])
		disp(['    \textbf{R044}    &     ',num2str(nUnits.R044.s1.n),' (',num2str(nUnits.R044.s1.nQ),')   &   ',num2str(nUnits.R044.s2.n),' (',num2str(nUnits.R044.s2.nQ),')  &  ',num2str(nUnits.R044.s3.n),' (',num2str(nUnits.R044.s3.nQ),')  &   ',num2str(nUnits.R044.s4.n),' (',num2str(nUnits.R044.s4.nQ),')  &   ',num2str(nUnits.R044.s5.n),' (',num2str(nUnits.R044.s5.nQ),')  &   ',num2str(nUnits.R044.s6.n),' (',num2str(nUnits.R044.s6.nQ),')  &\textbf{',num2str(nUnits.R044.total),' (',num2str(nUnits.R044.totalQ),')}\\ '])
		disp(['    \textbf{R050}    &     ',num2str(nUnits.R050.s1.n),' (',num2str(nUnits.R050.s1.nQ),')   &   ',num2str(nUnits.R050.s2.n),' (',num2str(nUnits.R050.s2.nQ),')  &  ',num2str(nUnits.R050.s3.n),' (',num2str(nUnits.R050.s3.nQ),')  &   ',num2str(nUnits.R050.s4.n),' (',num2str(nUnits.R050.s4.nQ),')  &   ',num2str(nUnits.R050.s5.n),' (',num2str(nUnits.R050.s5.nQ),')  &   ',num2str(nUnits.R050.s6.n),' (',num2str(nUnits.R050.s6.nQ),')  &\textbf{',num2str(nUnits.R050.total),' (',num2str(nUnits.R050.totalQ),')}\\ '])
		disp(['    \textbf{R064}    &     ',num2str(nUnits.R064.s1.n),' (',num2str(nUnits.R064.s1.nQ),')   &   ',num2str(nUnits.R064.s2.n),' (',num2str(nUnits.R064.s2.nQ),')  &  ',num2str(nUnits.R064.s3.n),' (',num2str(nUnits.R064.s3.nQ),')  &   ',num2str(nUnits.R064.s4.n),' (',num2str(nUnits.R064.s4.nQ),')  &   ',num2str(nUnits.R064.s5.n),' (',num2str(nUnits.R064.s5.nQ),')  &   ',num2str(nUnits.R064.s6.n),' (',num2str(nUnits.R064.s6.nQ),')  &\textbf{',num2str(nUnits.R064.total),' (',num2str(nUnits.R064.totalQ),')}\\ '])
        disp('    \hline')
        disp('    & & & & & & & &\\ % spacing hack')
        disp('  \end{tabular}')
        disp(' ')       
        disp('\label{tab:nUnits}')
        disp('\end{table}')
        
    otherwise
        error('Unrecognized tablemode. Better check that spelling.')   
end

disp(' '); 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('       END OF LATEX CODE'); 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'); disp(' ')

%% Finish up

% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

if cfg.writeFiles; diary off; end

cd(originalFolder)