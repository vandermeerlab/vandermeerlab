function [] = AMPX_spike_pipeline(cfg, varargin)
%% AMPX_spike_pipeline: This is a simple preprocessing pipeline for spike
%sorting the Naris experiment.  It will run the AMPX_to_ntt_JG script to
%create all the .ntt files for each of the four phases. It will then move
%these files to their own directories within the data dir, add a new FD
%file, and then it will run klustakwik for each session.  A log .txt file
%will be written.
%
% inputs:
% if only a specified phase (ie 'right' only) then this can be written as
% cfg.phases = {'right}.  Only the phases speficied here will be processed.


%% initialze some variables.
fname = cfg.fname;
cfg.phases = {'pre','left','right', 'post'};
cfg.tts_to_process = [5:8, 9:13, 15];
current_dir =(['G:\Naris\' fname(1:4) '\' fname(1:15)]);
cd(current_dir);
cfg.cd = current_dir;
threshold = -85;
cfg.thresh = [threshold threshold threshold threshold];
extract_varargin

% make the log file
fid = fopen(['Log_' datestr(now, 'dd_mm_yyyy_HH_MM') '.txt'], 'wt');

fprintf('\nPhases to be processed:\n')
for ii = 1:length(cfg.phases)
    fprintf([cfg.phases{ii} '\n'])
end
%% Create the ntt files for each session


for iphase = 1%:length(cfg.phases)
    fprintf(['Precprocessing ' cfg.phases{iphase} '...\n'])
    cfg.fname = [fname '-' cfg.phases{iphase} '.dat'];
    AMPX_to_ntt_EC(cfg)
    fprintf(fid, [datestr(now),cfg.phases{iphase} ' ntt Files created for channels' cfg.tts_to_process ' \n'])
    mkdir([cfg.phases{iphase} '_ntt\FD'])
    for ii = 1:length(cfg.tts_to_process)
        movefile([current_dir '\' fname '-' cfg.phases{iphase} '-TT' num2str(cfg.tts_to_process(ii)) '.ntt'] ,[current_dir '\' cfg.phases{iphase} '_ntt']);
        movefile([current_dir '\' fname '-' cfg.phases{iphase} '-TT' num2str(cfg.tts_to_process(ii)) '_pos.ntt'] ,[current_dir '\' cfg.phases{iphase} '_ntt']);

    end
    copyfile('C:\Users\mvdmlab\Dropbox\Matlab\GitHub\Si Analysis\Amplipex preprocesing workflow\Batch.txt',[current_dir '\' cfg.phases{iphase} '_ntt']);  % copy the batch file to the processing folder
    cd([current_dir '\' cfg.phases{iphase} '_ntt'])
    fprintf(fid, [datestr(now),cfg.phases{iphase} ' ntt folder complete\n'])
    
    RunClustBatch
    cd(current_dir)
    fprintf(fid, [datestr(now),cfg.phases{iphase} ' Klustakwik complete \n\n\n'])
    
end

end


