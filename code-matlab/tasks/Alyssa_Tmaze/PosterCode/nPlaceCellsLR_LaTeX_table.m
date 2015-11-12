%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%    Generate LaTeX TableSpeak for number of L&R-only Place Cells     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Tmaze data

% ABOUT
% This script loads the output from the script nUnitsRecorded_GenData and
% produces command window text that generates a LaTeX table. The text must
% be pasted into the LaTeX document manually.

% The user can choose to have the command window text saved to a text file
% (MATLAB'S diary function) so that each run is saved permanently.

% aacarey Sept 2015

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 WHAT DO YOU WANT THIS SCRIPT TO DO?                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];

% How to plot the table?
cfg.tablemode = 'default'; % future possibilities? 'poster' or 'paper' (not implemented)

% Which file to load and where to find it?
cfg.input_fn = 'nPlaceCellsLR'; % the script knows the proper file extension (.mat) so don't put it here
cfg.input_fd = 'D:\My_Documents\TmazePaper\data';

% Would you like to save the command window text (diary)?
cfg.writeDiary = 1; % 1 if yes, 0 if no

% Where to send the diary? (applies only if cfg.writeDiary = 1)
cfg.output_fd = 'D:\My_Documents\TmazePaper\visuals'; 

% Customize your diary name (applies only if cfg.writeDiary = 1)
cfg.prefix = ''; % prepend to filename
cfg.output_fn = 'nPlaceCellsLR_LaTeX'; % 'nUnitsRecorded_LaTeX'
cfg.suffix = ''; % append to filename. 

%% Load that stuff

iWasHere = pwd;

cd(cfg.input_fd)
load([cfg.input_fn,'.mat'])

if cfg.writeDiary
    diary on
    cd(cfg.output_fd)
    diary([cfg.output_fn,'.txt'])
    cd(cfg.output_fd)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(' ')
    disp(datestr(now)); disp(' ') % show me the date and the time i ran this
    disp(['   Loading ',cfg.input_fn]); disp(' ')
    disp(['   This file was created on ',nUnits.datetime]); disp(' ')
    disp(['   You have chosen tablemode [',cfg.tablemode,']']); disp(' ')
    disp(' ')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

%% I like it when you speak to me in LaTex
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('      LATEX CODE FOLLOWS:'); 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'); disp(' ')


switch cfg.tablemode
    case 'default'
disp('\begin{table}[h]')
disp('\centering')
disp('\caption[Neural units with fields on the track arms]{\textmd{Neural units with fields on the track arms for three rats across each of their six sessions. (L) left arm units; (R) right arm units}}')
disp('\small')
disp('	\begin{tabular}{ l c c c c c c c c c c c c c c c c c c c c}')
disp('		\hline		')
disp('		  & \multicolumn{2}{c}{\textbf{Session 1}}    &&   \multicolumn{2}{c}{\textbf{Session 2}}   & & \multicolumn{2}{c}{\textbf{Session 3}} && \multicolumn{2}{c}{\textbf{Session 4}}  && \multicolumn{2}{c}{\textbf{Session 5}}  && \multicolumn{2}{c}{\textbf{Session 6}}  && \multicolumn{2}{c}{\textbf{Total}} \\')
disp('		\hline')
disp('		\textbf{Rat ID} & \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} && \textbf{L} & \textbf{R} \\')
disp('		\hline')
disp(['		\textbf{R042}  &   ',num2str(nPC.R042.s1.L),'  &   ',num2str(nPC.R042.s1.R),'  &&  ',num2str(nPC.R042.s2.L),'  &  ',num2str(nPC.R042.s2.R),'  &&  ',num2str(nPC.R042.s3.L),'  &  ',num2str(nPC.R042.s3.R),'  &&  ',num2str(nPC.R042.s4.L),'  &   ',num2str(nPC.R042.s4.R),'  &&  ',num2str(nPC.R042.s5.L),'  &  ',num2str(nPC.R042.s5.R),'  &&  ',num2str(nPC.R042.s6.L),'  &  ',num2str(nPC.R042.s6.R),'  &&  ','\textbf{',num2str(nPC.R042.total.L),'}  &  ','\textbf{',num2str(nPC.R042.total.R),'}\\']) 
disp(['		\textbf{R044}  &   ',num2str(nPC.R044.s1.L),'  &   ',num2str(nPC.R044.s1.R),'  &&   ',num2str(nPC.R044.s2.L),'  &   ',num2str(nPC.R044.s2.R),'  &&  ',num2str(nPC.R044.s3.L),'  &  ',num2str(nPC.R044.s3.R),'  &&  ',num2str(nPC.R044.s4.L),'  &   ',num2str(nPC.R044.s4.R),'  &&   ',num2str(nPC.R044.s5.L),'  &   ',num2str(nPC.R044.s5.R),'  &&   ',num2str(nPC.R044.s6.L),'  &  ',num2str(nPC.R044.s6.R),'  &&  ','\textbf{',num2str(nPC.R044.total.L),'}  &   ','\textbf{',num2str(nPC.R044.total.R),'}\\']) 
disp(['		\textbf{R050}  &  ',num2str(nPC.R050.s1.L),'  &  ',num2str(nPC.R050.s1.R),'  &&  ',num2str(nPC.R050.s2.L),'  &  ',num2str(nPC.R050.s2.R),'  &&  ',num2str(nPC.R050.s3.L),'  &  ',num2str(nPC.R050.s3.R),'  &&  ',num2str(nPC.R050.s4.L),'  &  ',num2str(nPC.R050.s4.R),'  &&  ',num2str(nPC.R050.s5.L),'  &  ',num2str(nPC.R050.s5.R),'  &&  ',num2str(nPC.R050.s6.L),'  &  ',num2str(nPC.R050.s6.R),'  &&  ','\textbf{',num2str(nPC.R050.total.L),'}  &  ','\textbf{',num2str(nPC.R050.total.R),'}\\']) 
disp(['		\textbf{R064}  &  ',num2str(nPC.R064.s1.L),'  &  ',num2str(nPC.R064.s1.R),'  &&  ',num2str(nPC.R064.s2.L),'  &  ',num2str(nPC.R064.s2.R),'  &&  ',num2str(nPC.R064.s3.L),'  &  ',num2str(nPC.R064.s3.R),'  &&  ',num2str(nPC.R064.s4.L),'  &  ',num2str(nPC.R064.s4.R),'  &&  ',num2str(nPC.R064.s5.L),'  &  ',num2str(nPC.R064.s5.R),'  &&  ',num2str(nPC.R064.s6.L),'  &  ',num2str(nPC.R064.s6.R),'  &&  ','\textbf{',num2str(nPC.R064.total.L),'}  &  ','\textbf{',num2str(nPC.R064.total.R),'}\\']) 
disp('		\hline')
disp('		& & & & & & & &\\ % spacing hack')
disp('		& & & & & & & &\\ % spacing hack')
disp('	\end{tabular}')
disp(' ')
disp('	\label{tab:nPCells}')
disp('\end{table}')
        
    otherwise
        error('Unrecognized tablemode. Better check that spelling.')   
end

disp(' '); 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('       END OF LATEX CODE'); 
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'); disp(' ')

disp(['Total left cells from all rats: ',num2str(nPC.TotalAllRats.L)]);
disp(['Total right cells from all rats: ',num2str(nPC.TotalAllRats.R)]);
disp(['Total place cells from all rats: ',num2str(nPC.TotalAllRats.L + nPC.TotalAllRats.R)]);disp(' '); disp(' '); disp(' ')

%% Close the doors

if cfg.writeDiary
    diary off
end

cd(iWasHere)