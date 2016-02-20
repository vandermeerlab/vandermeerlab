%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%    Generate figure for thesis that has repr PSDs for candiates      %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% relies on fieldtrip functions
% plots power for candidate events and fake events (running) and saves image files
% does not save any variables generated
% loads files named *candidates.mat only

% ACarey May 2015

%%
clear

%% Choose settings
% master config for this script
cfg = [];
cfg.sessionFD   = {'E:\data\promoted\R042\R042-2013-08-16','E:\data\promoted\R044\R044-2013-12-20','E:\data\promoted\R050\R050-2014-04-03','E:\data\promoted\R064\R064-2015-04-23'}; % choose 3 representative sessions (one for each rat)
cfg.whichEvents = 'all'; %'prerecord', 'postrecord','task'
cfg.foi         = 50:300; %1:360; % the frequencies of interest (Hz)
cfg.toi         = -0.125:0.001:0.125;% -0.3:0.005:0.3; % the times of interest (before and after the center time for each event)

cfg.output_fd   = 'E:\Documents\TmazePaper\data';
cfg.output_fn   = 'TFRdata';

cfg.writeDiary = 0;
cfg.prefix      = ''; % prefix for image files
cfg.suffix      = ['-',date]; % suffix for image files; does NOT include extension

cfg.warningoff  = 1; % 1 don't show 8500 warnings that ft displays; else show warnings

%% check stuff
cfg_temp = []; 
cfg_temp.requireExpKeys = 1;
cfg_temp.ExpKeysFields = {'goodSWR','goodTheta'};
cfg_temp.requireMetadata = 1;
cfg_temp.MetadataFields = {'taskvars'};
cfg_temp.requireCandidates = 1;

proceed = checkTmazeReqs(cfg_temp);
cfg_temp = [];
rats = {'R042','R044','R050','R064'};

if proceed
    
    originalFolder = pwd;
    if cfg.warningoff
    warning off % because otherwise the avg colour of command window is orange
    end
    
    % do the thing
    for iFD = 1:length(cfg.sessionFD)
        cd(cfg.sessionFD{iFD})
        workingFolder = pwd;
        [~,sessionID,~] = fileparts(pwd);
        savename = [cfg.prefix,'ThreePSDs',cfg.suffix];
        if cfg.writeDiary % save command window text
            cd(cfg.output_fd)
            diary([savename,'.txt'])
            cd(cfg.sessionFD{iFD})
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
        
        % if there was HS detachment, we need to restrict csc
        [~,sessionID,~] = fileparts(pwd);
        if strcmp(sessionID(1:4),'R044')
            % 1. restrict based on HS_detach_times. 
            if strcmp(sessionID,'R044-2013-12-21') || strcmp(sessionID,'R044-2013-12-22')
                disp('Excluding HS detach times')
                load(FindFile('*HS_detach_times.mat'))
                temp_iv = iv(t_start/10^6,t_end/10^6);
                cscR = restrict(csc,temp_iv);%%%%%%%%%%%%%%%%%%
            end
            % 2. restrict based on metadata.detachIV
            disp('Restricting CSC')
            cscR = restrict(csc,metadata.detachIV);
        end
        
        % since there are gaps from starting/stopping the recording, need to
        % interpolate gaps and fill with NaNs or something:
        data = ft_read_neuralynx_interp(fc); % this annoyingly prints 100+ lines of warning if warning not off
        
        
        % choose which candidates to use
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
        
        % generate ft trials from the events
        
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
        if strcmp(sessionID(1:4),'R044')
            trl_running = ft_faketrlFromIV(cfg_temp,metadata.taskvars.trial_iv,cscR.tvec,data.hdr);
        else
            trl_running = ft_faketrlFromIV(cfg_temp,metadata.taskvars.trial_iv,csc.tvec,data.hdr);
        end
        cfg_temp = [];
        cfg_temp.trl = trl_running;
        data_trl_running = ft_redefinetrial(cfg_temp,data);
        
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
        cfg_tfr.toi          = cfg.toi; % times of interest % ex: 2s before:step_size:+2s after
        
        % if NaNs appear in TF: if your data is defined e.g. between -1 and +1 seconds, and using
        % a cfg.t_ftimwin of 0.5 seconds, freqanalysis_mtmconvol will give non-NaN output only
        % between -0.75 and 0.75 seconds

        TFR = ft_freqanalysis(cfg_tfr, data_trl2);
       
       
        %use same config as above ^^
        
        TFR_running = ft_freqanalysis(cfg_tfr, data_trl_running);
        
        TFRdata.(rats{iFD}).TFR = TFR;
        TFRdata.(rats{iFD}).TFR_running = TFR_running;
        
        if cfg.writeDiary, diary off, end
    end
    
    %% save data file
    
    cd(cfg.output_fd)
    save(cfg.output_fn,'TFRdata')
    cd(originalFolder)
  
    if cfg.warningoff, warning on, end
    
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('~~~      End of run       ~~~') 
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
end