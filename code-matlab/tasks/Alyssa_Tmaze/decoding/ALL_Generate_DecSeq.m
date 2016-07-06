% basic batch script to generate sequence detection based on decoded data

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
    
    cfg_decSeq = [];
    cfg_decSeq.minSeqLength = 8;
    cfg_decSeq.removeInterneurons = 1;
    cfg_decSeq.nMinNeurons = 1;
    cfg_decSeq.output_file_prefix = 'R2_';
    
    for iFD = 1:length(fd)
        
        cd(fd{iFD});
        
        Generate_DecSeq(cfg_decSeq);
        
    end
end
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')