%%
cfg = []; cfg.rats = {'R064'};
cfg.requireMetadata = 1;
cfg.requireCandidates = 1;
fd = sort(getTmazeDataPath(cfg));

%% 
profile on
for iFD = 1:length(fd)
    
    cd(fd{iFD});
    mkdir('files')
    
%     cfg = [];
%     cfg.output_file_prefix = 'CS_PRE_'; % prefix use when writing files
%     cfg.whichEvents = 'prerecord';
%     
%     Generate_CorrScores(cfg);
       
    cfg .output_file_prefix = 'YTPOSTER_'; % prefix use when writing files
    cfg.whichEvents = 'all';
    cfg.matchFields = 1;
    
    Generate_CorrScores(cfg);
    
end

profile viewer