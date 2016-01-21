%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%         COUNT PLACE CELLS USED IN COOCCURRENCE ANALYSIS             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tmaze project
%
% This script counts how many place cells were used in CoOccurrence
% analysis. It can also be used to count how many cells we have of a
% certain category, such as those with ratings of 3 or higher that have
% place fields anywhere along the left trajectory.
% LaTeX code for a table is printed to the command window
%
% aacarey Jan 2016 (from thesis script nPlaceCells)
%
clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% Location of inputData
cfg.input_fd = [pwd,'\data'];

% Do you want to save the command window text?
cfg.writeFiles = 1;

% Where to save the text file?
cfg.output_fd = cfg.input_fd;

% Which rats?
cfg.rats = TmazeRats;

% Which S was used for co-occurrence analysis?
cfg.whichS = 'unique'; %'unique' are cells with place fields on only one arm; 
% 'nonunique' are cells that were active on a particular arm but may have 
% fields on the opposite arm; 'trajectory' are cells with fields along the 
% L or R trajectory

% Select which cells to keep (these are direct inputs into SelectTS targeting S.usr.rating)
cfg.operation = '<=';
cfg.threshold = 5;

%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.temp);
end

%% Info and data loading section

originalFolder = pwd;

if cfg.writeFiles % save command window text
    cd(cfg.output_fd)
    diary(['nPlaceCells_text','.txt'])
    cd(cfg.output_fd)
    disp(' ')
end
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                   COUNTING PLACE CELLS                      ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

fprintf('\nDATE of script run: %s\n',datestr(now));

cd(cfg.input_fd)

load('inputData.mat')
fprintf('\nThis inputData was generated on %s\n',inputData.hdr.date);

fprintf('\nCounting %s place cells with ratings %s %d \n',cfg.whichS,cfg.operation,cfg.threshold);

disp(' ')
switch cfg.whichS
    case 'unique'
        cfg.whichS = 'S_arm_unique';
    case 'nonunique'
        cfg.whichS = 'S_arm';
    case 'trajectory'
        cfg.whichS = 'S_traj';
    otherwise
        diary off
        error('Unrecognized cfg.whichS. Better check that spelling ^_^')
end

arms = {'L','R'};

%% Generate LaTex code

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('      LATEX CODE FOLLOWS:'); 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'); disp(' ')

disp('\begin{table}[h]')
disp('\centering')
disp('\caption[Neural units with fields on the track arms]{\textmd{Neural units with fields on the track arms for three rats across each of their six sessions. (L) left arm units; (R) right arm units}}')
disp('\small')
disp('	\begin{tabular}{ l c c c c c c c c c c c c c c c c c c c}')
disp('		\hline		')
disp('		  & \multicolumn{2}{c}{\textbf{Session 1}}    &&   \multicolumn{2}{c}{\textbf{Session 2}}   & & \multicolumn{2}{c}{\textbf{Session 3}} && \multicolumn{2}{c}{\textbf{Session 4}}  && \multicolumn{2}{c}{\textbf{Session 5}}  & &\multicolumn{2}{c}{\textbf{Session 6}}  && \textbf{Total} \\')
disp('		\hline')
disp('		\textbf{Rat ID} & \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{LR} \\')
disp('		\hline')

% config for SelectTS
cfg_temp.threshold = cfg.threshold; cfg_temp.operation = cfg.operation; cfg_temp.verbose = 0;

% move through inputData and print latex code to command window
for iRat = 1:length(cfg.rats)
    
    nTotal = 0; % for collecting total place cells per rat
    fprintf('\t\t\t\\textbf{%s}',cfg.rats{iRat});
    
    for iSession = 1:length(inputData.(cfg.rats{iRat}))
        
        for iArm = 1:length(arms)
            
            S = inputData.(cfg.rats{iRat})(iSession).(arms{iArm}).(cfg.whichS);
            
            if cfg.threshold < 5 % then eliminate some cells
                S = SelectTS(cfg_temp,S,'rating');
            end
            
            nTotal = nTotal + length(S.t);
            fprintf('   &   %d',length(S.t))
            
        end % of arms
        
        fprintf('   &&')
        
    end % of sessions
    
    disp(['   \textbf{',num2str(nTotal),'}\\']) 
    
end % of rats

disp('		\hline')
disp('		& & & & & & & &\\ % spacing, this might not be necessary')
disp('		& & & & & & & &\\ % spacing')
disp('	\end{tabular}')

disp('	\label{tab:nPCells}')
disp('\end{table}')
    
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