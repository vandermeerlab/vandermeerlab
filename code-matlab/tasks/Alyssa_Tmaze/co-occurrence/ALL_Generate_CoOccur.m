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
proceed = checkTmazeReqs(cfg); % make sure we have everything

if proceed
    fd = sort(getTmazeDataPath(cfg));

    for iFD = 1:length(fd)
        
        cd(fd{iFD});
        
        cfg = [];
        cfg.output_file_prefix = 'CC_PRE_'; % prefix use when writing files
        cfg.whichEvents = 'prerecord';
        cfg.writeFiles = 1;
        cfg.plotOutput = 0;
        cfg.load_questionable_cells = 1;
        
        Generate_CoOccur(cfg);
        
    end
end