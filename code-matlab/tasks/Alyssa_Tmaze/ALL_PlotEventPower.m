%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                        Plot event power                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% relies on fieldtrip functions
% plots power for candidate events and fake events (running) and saves image files
% does not save any variables generated
% loads files named *candidates.mat only

% ACarey Mar 2015, with pieces from MvdM code 
%  -- edit May 2015

%%
clear

%% Choose settings
% master config for this script
cfg = [];
cfg.rats = {'R042','R044','R050'};
cfg.whichEvents = 'all'; %'prerecord', 'postrecord','task'
cfg.foi = 1:250; % the frequencies of interest (Hz)
cfg.writeFiles = 1;
cfg.output_fd = 'files'; % file directory for saving image files
cfg.prefix = ''; % prefix for image files
cfg.suffix = ''; % suffix for image files; does NOT include extension
cfg.writeDiary = 1; % save command window text 
cfg.warningoff = 1; % don't show 8500 warnings that ft displays; else show warnings

%% check stuff
cfg_temp = []; 
cfg_temp.requireExpKeys = 1;
cfg_temp.ExpKeysFields = {'goodSWR','goodTheta'};
cfg_temp.requireMetadata = 1;
cfg_temp.MetadataFields = {'taskvars'};
cfg_temp.requireCandidates = 1;
if cfg.writeDiary || cfg.writeFiles, cfg_temp.requireFiles = 1; end

proceed = checkTmazeReqs(cfg_temp);
cfg_temp = [];

if proceed
    
    cfg_temp.rats = cfg.rats;
    fd = sort(getTmazeDataPath(cfg_temp)); % get all session directories
    
    originalFolder = pwd;
    if cfg.warningoff
    warning off % because otherwise the avg colour of command window is orange
    end
    
    % do the thing
    for iFD = 1:length(fd)
        cd(fd{iFD})
        workingFolder = pwd;
        [~,sessionID,~] = fileparts(pwd);
        savename = [cfg.prefix,sessionID,'-EventPower',cfg.suffix];
        if cfg.writeDiary % save command window text
            cd([fd{iFD},'\files'])
            diary([savename,'.txt'])
            cd(fd{iFD})
            disp(date)
        end
        cprintf(-[0 0 1],['Working on',[' ',sessionID]]); disp(' ')
        
        disp('Loading ExpKeys')
        LoadExpKeys
        
        disp('Loading metadata')
        LoadMetadata
        %metaload = [];
        %metaload.suffix = '';
        %LoadMetadata3(metaload) % to get trial intervals
        %clear metaload
        
        disp('Loading candidates')
        LoadCandidates
        
        cfg_temp = [];
        fc = ExpKeys.goodSWR(1); % for ft_read_neurablah...
        cfg_temp.fc = fc; % for LoadCSC to get tvec for ft_maketrlFromIV()
        csc = LoadCSC(cfg_temp);
        
        %% if there was HS detachment, we need to replace HS detach times with zeros in CSC
        [~,sessionID,~] = fileparts(pwd);
        if strcmp(sessionID(1:4),'R044')
            % 1. add zeros within HS_detach_times. 
            if strcmp(sessionID,'R044-2013-12-21') || strcmp(sessionID,'R044-2013-12-22')
                
                load(FindFile('*HS_detach_times.mat'))
                temp_iv = iv(t_start/10^6,t_end/10^6);
                insertHere = zeros(size(csc.tvec));
                for iInterval = 1:length(temp_iv.tstart)
                   ugh = csc.tvec > temp_iv.tstart(iInterval) & csc.tvec < temp_iv.tend(iInterval);
                   insertHere = insertHere + ugh;
                end
                   insertHere = ~logical(insertHere);
                   csc.data(insertHere) = 0;
            end
            % 2. add zeros within metadata.detachIV
            
            insertHere = zeros(size(csc.tvec));
            for iInterval = 1:length(metadata.detachIV.tstart)
                ugh = csc.tvec > metadata.detachIV.tstart(iInterval) & csc.tvec < metadata.detachIV.tend(iInterval);
                insertHere = insertHere + ugh;
            end
            insertHere = ~logical(insertHere);
            csc.data(insertHere) = 0;
        end
        
        % since there are gaps from starting/stopping the recording, need to
        % interpolate gaps and fill with NaNs or something:
        data = ft_read_neuralynx_interp(fc); % this annoyingly prints 100+ lines of warning if warning not off
        
        
        %% choose which candidates to use
        fprintf('nEvents before cfg.whichEvents restrict: %d\n',length(evt.tstart));
        
        switch cfg.whichEvents
            case 'all'
                %evt = evt; 
            case 'prerecord'
                evt = restrict(evt,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
            case 'task'
                evt = restrict(evt,ExpKeys.task(1),ExpKeys.task(2));
            case ' postrecord'
                evt = restrict(evt,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
        end
        fprintf('nEvents after cfg.whichEvents restrict: %d\n',length(evt.tstart));
        
        %% generate ft trials from the events
        
        twin = [-1 1];
        
        disp('Generating trials from candidates')
        cfg_temp = [];
        cfg_temp.twin = twin;
        trl = ft_maketrlFromIV(cfg_temp,evt,csc.tvec,data.hdr);
              
        cfg_temp = [];
        cfg_temp.trl = trl;
        data_trl2 = ft_redefinetrial(cfg_temp,data);
        
        disp('Generating trials from fake events')
        cfg_temp = [];
        cfg_temp.num = length(evt.tstart);
        cfg_temp.twin = twin;
        trl_running = ft_faketrlFromIV(cfg_temp,metadata.taskvars.trial_iv,csc.tvec,data.hdr);
        
        cfg_temp = [];
        cfg_temp.trl = trl_running;
        data_trl_running = ft_redefinetrial(cfg_temp,data);
        
        %% Generate plot for events
        
        disp('Performing frequency analysis')
        
        cfg_tfr              = []; % start with empty cfg
        cfg_tfr.output       = 'pow';
        cfg_tfr.channel      = data_trl2.label;
        cfg_tfr.method       = 'mtmconvol';
        cfg_tfr.taper        = 'hanning';
        cfg_tfr.foi          = cfg.foi; %1:50; %1:250; %7:0.5:10; % frequencies of interest
        %cfg.t_ftimwin    = ones(size(cfg.foi)).*0.5; % window size: fixed at 0.5s
        cfg_tfr.t_ftimwin    = 5./cfg_tfr.foi;
        %cfg.tapsmofrq = cfg.foi*0.4;
        cfg_tfr.toi          = -0.5:0.005:0.5; % times of interest % ex: 2s before:step_size:+2s after
        
        % if NaNs appear in TF: if your data is defined e.g. between -1 and +1 seconds, and using
        % a cfg.t_ftimwin of 0.5 seconds, freqanalysis_mtmconvol will give non-NaN output only
        % between -0.75 and 0.75 seconds

        TFR = ft_freqanalysis(cfg_tfr, data_trl2);
        
        disp(' ')
        disp('Plotting power for candidate events')
        fh = figure;
        cfg_temp = []; cfg_temp.channel = data_trl2.label;
        ft_singleplotTFR(cfg_temp, TFR);
        grid on
        ylabel('Frequency (Hz)')
        if cfg.writeFiles && ~isempty(cfg.output_fd)
            cd([pwd,'\',cfg.output_fd])
            filename = [cfg.prefix,sessionID,'-CandEventPower',cfg.suffix,'.png'];
            set(fh, 'InvertHardCopy', 'off');
            print(fh,'-r300','-dpng',filename);
            close(fh);
            cd(workingFolder)
        end
       
        %% generate plot for running

        %use same config as above ^^
       
        TFR_running = ft_freqanalysis(cfg_tfr, data_trl_running);
        
        disp(' ')
        disp('Plotting power for fake events')
        figure
        cfg_temp = []; cfg_temp.channel = data_trl_running.label;
        ft_singleplotTFR(cfg_temp, TFR_running);
        grid on
        ylabel('Frequency (Hz)')
        if cfg.writeDiary, diary off, end
        if cfg.writeFiles && ~isempty(cfg.output_fd)
            cd([pwd,'\',cfg.output_fd])
            filename = [cfg.prefix,sessionID,'-FakeEventPower',cfg.suffix,'.png'];
            set(fh, 'InvertHardCopy', 'off');
            print(fh,'-r300','-dpng',filename);
            close(fh);
            cd(workingFolder)
        end
    end
    %%
    if cfg.warningoff, warning on, end
    
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~      End of run       ~~~') 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
end