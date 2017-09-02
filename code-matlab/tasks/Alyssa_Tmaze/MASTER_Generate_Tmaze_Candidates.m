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

% ACarey May 2015

%%

prefix = ''; % ex: "test_', will be prepended to the filename
suffix = '_HT1-3_noTheta'; %[date,'']; appended to filename

writeDiary = 1; % if 1, save command window text to session 'files' folders;
%                  the command window text includes things like "n events
%                  left after theta thresholding" and so on.

gen = [];
gen.load_questionable_cells = 0;
gen.SWRmethod = 'HT'; 
gen.MUAmethod = 'none';
gen.ThetaThreshold = [];

%% verify requisites: detection is a long process, don't want it erroring partway through 

cfg = []; 
cfg.requireExpKeys = 1;
cfg.ExpKeysFields = {'goodSWR','goodTheta'};
cfg.requireMetadata = 1;
cfg.MetadataFields = {'SWRfreqs'};
cfg.requireVT = 1;
cfg.requireHSdetach = 1;
if writeDiary, cfg.requireFiles = 1; end

proceed = checkTmazeReqs(cfg); 

%%
cfg.rats = {'R042','R044','R050','R064'};
originalFolder = pwd;

if proceed % Detect candidate replay events and save them as a .mat file

    fd = sort(getTmazeDataPath(cfg)); % get all session directories  
   
    for iFD = 1:length(fd)
        cd(fd{iFD});
        
        [~,session,~] = fileparts(fd{iFD});
        savename = strcat(prefix,session,'-candidates',suffix); % string concatenation
        
        if writeDiary % save command window text
            cd([fd{iFD},'\files'])
            diary([savename,'-final','.txt'])
            cd(fd{iFD})
            disp(date)
        end
        
        % tell me what you're working on
        cprintf(-[0 0 1],['Working on ',session]); % disp(['Working on ',session])
        disp(' ');
    
        % generate candidates 
        evt = GenCandidateEvents(gen); % T-maze specific
        
        if writeDiary, diary off, end
        
        % save candidates
        save([savename,'.mat'],'evt'); 
    end
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~ End of candidates run ~~~') 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
end

cd(originalFolder)
