%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%              COUNT NUMBER OF UNITS WE RECORDED                      %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Tmaze data

% ABOUT
% This script loads .t and ._t files and counts all the units recorded per
% session for each rat. The data is send to a directory specified by the
% user.
% It prints text to the command window and saves it in a text file, using
% MATLAB's diary function, for quick reference if needed. The text file is 
% saved in the same directory as the data.

% To generate a table in LaTeX speak, run the script nUnitsRecorded_LaTeX_table

% IMPORTANT
% Because of the script author's coding limitations, each rat must have 
% exactly 6 session folders.

% aacarey Sept 2015

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

cfg = [];

% Where to send the data and diary?
cfg.output_fd = 'D:\My_Documents\TmazePaper\data'; 

% Which rats to include?
cfg.rats = {'R042','R044','R050','R064'};

% Should I load the 5's?
cfg.load_questionable_cells = 1;

% Would you like to save the data and the diary?
cfg.writeData = 1; 
cfg.writeDiary = 1; 

% Customize your filenames (same for data and diary)
cfg.prefix = ''; % prepend to filename
cfg.output_fn = 'nUnitsRecorded';
cfg.suffix = ''; % append to filename. 

%%
if cfg.writeDiary % save command window text
    cd(cfg.output_fd)
    diary([cfg.prefix,cfg.output_fn,cfg.suffix,'.txt'])
    cd(cfg.output_fd)
end

% set up LoadSpikes' config
please = [];
please.load_questionable_cells = cfg.load_questionable_cells;
please.useClusterFile = 0;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp(datestr(now)); disp(' ') % show me the date and the time i ran this
disp('   The following data include these rats: '); disp(cfg.rats)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

originalFolder = pwd; % remember where we started

% Make something to collect each rat's session counts into. s1 is session 1
% and it will hold the number of units from that session.
empty = struct('s1',[],'s2',[],'s3',[],'s4',[],'s5',[],'s6',[],'total',[]);

TotalAllRats = 0; % grand number of units for experiment

% go through each rat and collect the data
for iRat = 1:length(cfg.rats);
    
    collect = empty; % make a copy of empty to put stuff into
    total = 0; % this will collect the number of units in total for a single rat

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
        
        % collect total cells recorded for all rats
        TotalAllRats = TotalAllRats + n;
        
        % collect total for each rat     
        total = total + n;
        
        % collect cells from each session
        collect.(['s',num2str(iFD)]).n = n; %  wow eh. 
        % ^ tells matlab to assign the value to (for ex), empty.s1 for the first file directory
        % ['s',num2str(iFD)] if run, generates s1 or s2 and so on
        collect.(['s',num2str(iFD)]).sessionID = sessionID; % if needed for future data handling
        
        
    end
    collect.total = total; % the total cells recorded for the current rat
    nUnits.(cfg.rats{iRat}) = collect; % makes, for example, nUnits.R042 = collect;
    disp(' ')
    disp(['***Total units for ',cfg.rats{iRat},': ',num2str(total)]); disp(' ')
end

nUnits.TotalAllRats = TotalAllRats;
nUnits.datetime = datestr(now);

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp(['  Total cells for all cfg.rats: ',num2str(TotalAllRats)])
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

%% Tie loose ends

if cfg.writeDiary; diary off; end

if cfg.writeData
    cd(cfg.output_fd)
    save([cfg.prefix,cfg.output_fn,cfg.suffix],'nUnits')
end

cd(originalFolder)