% basic batch script to generate co-occurrence data

%%
cfg = []; 
cfg.requireExpKeys = 1;
cfg.ExpKeysFields = {'prerecord','postrecord','goodTheta','pathlength'};
cfg.requireMetadata = 1;
cfg.MetadataFields = {'coord','taskvars'};
cfg.requireCandidates = 1;
cfg.requireVT = 1;
cfg.requireTimes = 1; %R042 only
cfg.requireFiles = 1;
cfg.rats = {'R042','R044','R050','R064'};
proceed = checkTmazeReqs(cfg); % make sure we have everything

if proceed
    fd = sort(getTmazeDataPath(cfg));
    cfg = [];
    cfg.whichEvents = 'all';
    
    for iFD = 1:length(fd)
        
        cd(fd{iFD});
        
        
        switch cfg.whichEvents
            case 'all'
                cfg.output_file_prefix = 'CC_ALL_';
            case 'prerecord'
                cfg.output_file_prefix = 'CC_PRE_';
            case 'postrecord'
                cfg.output_file_prefix = 'CC_POST_';
            case 'taskrest'
                cfg.output_file_prefix = 'CC_TASKREST_';
            case 'task'
                cfg.output_file_prefix = 'CC_TASK_';
            otherwise
                error('Unrecognized cfg.whichEvents')
        end
                
        cfg.writeFiles = 1;
        cfg.plotOutput = 0;
        cfg.load_questionable_cells = 1;
        cfg.outputPConly = 1;
        
        Generate_CoOccur(cfg);
        
    end
end
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')