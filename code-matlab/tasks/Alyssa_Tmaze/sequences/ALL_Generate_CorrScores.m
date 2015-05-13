%%
cfg = []; cfg.rats = {'R042','R050'};
cfg.requireMetadata = 1;
cfg.requireCandidates = 1;
fd = sort(getTmazeDataPath(cfg));

%% 
tic
for iFD = 12%:length(fd)
    
    cd(fd{iFD});
    
%     cfg = [];
%     cfg.output_file_prefix = 'CS_PRE_'; % prefix use when writing files
%     cfg.whichEvents = 'prerecord';
%     
%     Generate_CorrScores(cfg);
    
    cfg.output_file_prefix = 'CS_M_ALL_'; % prefix use when writing files
    cfg.whichEvents = 'all';
    cfg.matchFields = 1;
    
    Generate_CorrScores(cfg);
    
end

toc