%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                   %%%%
%%%%             MASTER Generate Tmaze Candidate Events                %%%%
%%%%                                                                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the pool of candidates that will be used for further
% analysis; relies on the function GenCandidateEvents()

% SWR & MUA detection
% Speed & Theta thresholding
% NActive cells thresholding
% ResizeIV (interval expansion/contraction) optional

% aacarey May 2015, initial version
% aacarey 02 MAR 2019, added extra input config options and coordination
%         with main analysis script. Can still be run independently.

clearvars -except CFG
%% CHOOSE WHAT YOU WANT THIS SCRIPT TO DO

cfg.suffix = '_trymatch'; %[date,'']; appended to filename

cfg. writeDiary = 1; % if 1, save command window text to session 'files' folders;
%                  the command window text includes things like "n events
%                  left after theta thresholding" and so on.

cfg.SWRmethod = 'AM'; % 'AM' for amSWR (frequency content similarity), 'HT' for OldWizard (hilbert transform), 'TR' for photonic (transient detection), or 'none' to skip
cfg.MUAmethod = 'AM'; % 'AM' for amMUA, or 'none' to skip MUA detection
cfg.weightby = 'amplitude'; % this applies to 'AM' amSWR and 'TR' photonic, but not 'HT' OldWizard
cfg.load_questionable_cells = 1;
cfg.DetectorThreshold = 4; % the threshold you want for generating IV data edges
cfg.ThreshMethod = 'raw';
cfg.ThetaThreshold = 2;

%% Handle inputs from main analysis script, if exist

if exist('CFG','var') && ~isempty(CFG.output_folder)
   cfg.writeDiary = 1;
   CFG.candidates_filename_suffix = cfg.suffix; % for loading the correct candidates in later scripts
end

%% verify requisites: detection is a long process, don't want it erroring partway through 

cfg_check.requireExpKeys = 1;
cfg_check.ExpKeysFields = {'goodSWR','goodTheta'};
cfg_check.requireMetadata = 1;
cfg_check.MetadataFields = {'SWRfreqs'};
cfg_check.requireVT = 1;
cfg_check.requireHSdetach = 1;
if cfg.writeDiary, cfg.requireFiles = 1; end

proceed = checkTmazeReqs(cfg_check); 

%%
cfg.prefix = ''; % ex: "test_', will be prepended to the filename
cfg_path.rats = {'R042','R044','R050','R064'};
originalFolder = pwd;

if proceed % Detect candidate replay events and save them as a .mat file

    fd = sort(getTmazeDataPath(cfg_path)); % get all session directories  
   
    for iFD = 1:length(fd)
        cd(fd{iFD});
        
        [~,session,~] = fileparts(fd{iFD});
        savename = strcat(cfg.prefix,session,'-candidates',cfg.suffix); % string concatenation
        
        if cfg.writeDiary % save command window text
            cd([fd{iFD},'\files'])
            diary([savename,'-final','.txt'])
            cd(fd{iFD})
            disp(date)
        end
        
        % tell me what you're working on
        cprintf(-[0 0 1],['Working on ',session]); % disp(['Working on ',session])
        disp(' ');
    
        % generate candidates 
        cfg_gen.load_questionable_cells = cfg.load_questionable_cells;
        cfg_gen.SWRmethod = cfg.SWRmethod; 
        cfg_gen.MUAmethod = cfg.MUAmethod;
        cfg_gen.ThetaThreshold = cfg.ThetaThreshold;
        evt = GenCandidateEvents(cfg_gen); % T-maze specific
        
        if cfg.writeDiary, diary off, end
        
        % save candidates
        save([savename,'.mat'],'evt'); 
    end
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~ End of candidates run ~~~') 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
end

cd(originalFolder)
