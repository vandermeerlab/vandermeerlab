% basic batch script to generate co-occurrence data

%%
cfg = []; cfg.rats = {'R042','R044','R050'};
cfg.requireMetadata = 1;
cfg.requireCandidates = 1;
fd = sort(getTmazeDataPath(cfg));

%% 
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