function [csc_tsd] = load_CSC_spike_gadgets(cfg_in)
        % function csc_tsd = load_CSC_spike_gadgets(cfg_in)
        %
        % loads SpikeGadgets ntxchx.dat files which contain raw data for 1
        % channel
    
        % INPUTS: only cfg
        %
        % cfg.fc: cell array containing filenames to load
        %     if no file_list field is specified, loads all *.Ncs files in

        % TODO deal with the different filename for unique timestamp file
        % and possibly multiple lfp data files
        % TODO write warning that the names should be in the same folder
        % TODO check the uV and put a warning
        % TODO check that time and data are the same length 
    
        cfg_def.fc = {};
        cfg_def.fc = {'2023-03-22_r206_sq21.LFP_nt27ch3.dat'};
    %     cfg_def.TimeConvFactor = 1/30000; % 30 kHz sampling rate, maybe better to get it from the file?
        cfg_def.VoltageConvFactor = 1; % 1 means output in volts, 1000 in mV, 10^6 in uV
        cfg_def.decimateByFactor = [];
        cfg_def.verbose = 1;
        cfg_def.lfp_times_fn = {'2023-03-22_r206_sq21.timestamps.dat'};
    
        if ~exist('cfg_in', 'var') % temp
            cd('C:\shared\DATA\Screening_DM\r206\2023-03-22_r206_sq21\2023-03-22_r206_sq21.LFP\')
            cfg_in = cfg_def;
        end
        mfun = mfilename;
    
        cfg = ProcessConfig(cfg_def, cfg_in, mfun); % this takes fields 
        % from cfg_in and puts them into cfg
    
        if isempty(cfg.fc) % no filelist provided, load everything
            cfg.fc = FindFiles('*.dat');
        else    
            if ~isa(cfg.fc,'cell')
                error('LoadCSC: cfg.fc should be a cell array.');
            end
        end
    
        if isempty(FindFiles(cfg.fc{1}))
            error('File does not exist in directory. Check spelling.');
        end
    
        fc = sort(cfg.fc);
        nFiles = length(fc);
    
        % track sample counts for different files, throw error if not equal
        sample_count_tvec = nan(nFiles, 1);
        sample_count_data = nan(nFiles, 1);
    
        csc_tsd = tsd;
    
        % TODO check this according to what is in the header of the file
        if cfg.VoltageConvFactor == 1
            csc_tsd.units = 'V';
        elseif cfg.VoltageConvFactor == 1000
            csc_tsd.units = 'mV';
        elseif cfg.VoltageConvFactor == 10^6
            csc_tsd.units = 'uV';
        elseif cfg.VoltageConvFactor == 0
            csc_tsd.units = 'ADbits';
        else
            error('Input voltage conversion factor is not a valid value.')
        end
    
        if cfg.verbose
            fprintf('%s: Loading %d file(s)...\n', mfun, nFiles)
        end
    
        % 1. read the timestamps file (one file for all LFP files)
        disp('Reading LFP timestamps file...')
        these_LFP_times = readTrodesExtractedDataFile(cfg.lfp_times_fn);
        orig_Fs = these_LFP_times.clockrate;
    
        tvec = double(these_LFP_times.fields.data) / orig_Fs; % now in s
    
        % decimate data if specified
        if ~isempty(cfg.decimateByFactor)
            
            fprintf('%s: Decimating tvec by factor %d...\n', mfun, ...
                cfg.decimateByFactor)
            tvec = tvec(1:cfg.decimateByFactor:end);
            % TODO deal with the 'header' at some point
    %         hdr.SamplingFrequency = hdr.SamplingFrequency./cfg.decimateByFactor;
        end
    
        for iF = 1:nFiles
            fname = fc{iF};
    
            disp(['Reading LFP data file: ' fname '...'])
            % NOTE: the data should be in uV so to match the desired unit, 
            % we need to set it to V [TO IMPROVE ED]
            this_LFP =  readTrodesExtractedDataFile(fname);
            lfp_data = double(this_LFP.fields.data) * this_LFP.voltage_scaling;
            lfp_data = lfp_data/10^6;
            

            if isempty(lfp_data) % TODO check if this is ever the case
                warning(['No csc data (disabled tetrode channel). ' ...
                    'Consider deleting ',fname,'.']);
            end
    
            % track sizes
            sample_count_tvec(iF) = length(tvec);
            sample_count_data(iF) = length(lfp_data);
    
            % check if the data is the same length for each channel.  
            if iF >1 && length(lfp_data) ~= length(csc_tsd.data(iF-1,:))
                message = 'Data lengths differ across channels.';
                error(message);
            end

            %TODO do test of length of timestamp and lfp data should be the
            %same
    
            % decimate data if specified
            if ~isempty(cfg.decimateByFactor)
                
                fprintf('%s: Decimating by factor %d...\n', mfun, ...
                    cfg.decimateByFactor)
                lfp_data = decimate(lfp_data, cfg.decimateByFactor);
                tvec = tvec(1:cfg.decimateByFactor:end);
                % TODO deal with this header thing
    %             hdr.SamplingFrequency = hdr.SamplingFrequency./cfg.decimateByFactor;   
            end
    
            % done, add to tsd
            csc_tsd.tvec = tvec;
            csc_tsd.data(iF,:) = lfp_data;
            
            [~, fn, fe] = fileparts(fname);
            csc_tsd.label{iF} = cat(2, fn, fe);
    
            % TODO deal with the header
    %         csc_tsd.cfg.hdr{iF} = hdr;
    
        end
    % check if anything unequal --> error
    if numel(unique(sample_count_data)) > 1
        error('Data sizes unequal.');
    end
    
    if numel(unique(sample_count_tvec)) > 1
       error('tvec sizes unequal.');
    end
    
    % check if ExpKeys available
    keys_f = FindFiles('*keys.m');
    if ~isempty(keys_f)
        run(keys_f{1});
        csc_tsd.cfg.ExpKeys = ExpKeys;
    end
    
    % add sessionID
    [~, csc_tsd.cfg.SessionID, ~] = fileparts(pwd);
    
    % housekeeping
    csc_tsd.cfg.history.mfun = cat(1,csc_tsd.cfg.history.mfun,mfun);
    csc_tsd.cfg.history.cfg = cat(1,csc_tsd.cfg.history.cfg,{cfg});
end

%    
% % 
% %         % load raw LFP data for each file
% %         [Timestamps, ~, SampleFrequencies, NumberOfValidSamples, ...
% %             Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);
%          
% 
% 
%     end
% 
%     keyboard
% 
%     % initialisation
%     csc_tsd = [];
%     
%     % checks
%     if ~isfile([lfp_folder lfp_times_fn])
%         warning(['LFP times file not found at ' lfp_folder lfp_fn ...
%             '- might need to extract it?'])
%         return
%     end
% 
%     % read the data
% %     disp('Reading LFP timestamps file...')
% %     these_LFP_times = readTrodesExtractedDataFile([lfp_folder lfp_times_fn]);
% %     orig_Fs = these_LFP_times.clockrate;
% %     if isempty(lfp_new_Fs)
% %         Fs = orig_Fs;
% %     else
% %         Fs = lfp_new_Fs; % the sampling frequency that the file has been 
% %         % extracted to
% %     end
% %     lfp_t = double(these_LFP_times.fields.data) / orig_Fs; % now in s
% 
%     disp(['Reading LFP data file: ' lfp_fn '...'])
%     this_LFP =  readTrodesExtractedDataFile([lfp_folder lfp_fn]);
%     lfp_data = double(this_LFP.fields.data) * this_LFP.voltage_scaling;
%     % should be in... uvolts?
%     lfp_unit = 'uV';
% 
%     % more checks
%     if size(lfp_data,1) ~= size(lfp_t,1)
%         warning('LFP data and LFP times different??')
%         keyboard
%     end
% 
%     % if requested, downsample
%     if ds>0
%         [lfp_data, lfp_t] = resample(lfp_data, lfp_t, ds,'pchip'); % use interpolation and an anti-aliasing filter to resample the signal at a uniform sample rate
%         Fs = ds; % the new sampling rate
%         disp(['Downsampled lfp data to ' num2str(ds) 'Hz']);
%     end
% 
%     % convert LFP into mvdmlab TSD structure
%     % Prepare the tsd structure
%     cfg_in.fc = {lfp_fn};
%     mfun = mfilename; % get name of the current code to store it in the 
%     % history field of tsd
%     cfg_def.fc = {};
%     cfg_def.verbose = 1;
%     cfg = ProcessConfig(cfg_def, cfg_in, mfun);
% 
%     csc_tsd = tsd;
%     csc_tsd.data = lfp_data';
%         
%     csc_tsd.cfg.Fs = Fs; % I think this is sampling rate
%     csc_tsd.cfg.hdr{1}.Fs = Fs;
% 
%     % construct tvec & add to tsd
%     
%     csc_tsd.tvec = lfp_t';
% 
%     [~,fn,fe] = fileparts(lfp_fn);
%     csc_tsd.label =  cat(2, fn, fe);
%     csc_tsd.units = lfp_unit;
%     csc_tsd.reference = this_LFP.reference;
%     csc_tsd.cfg.SessionID = lfp_fn;
% 
%     % housekeeping
%     csc_tsd.cfg.history.mfun = cat(1,csc_tsd.cfg.history.mfun,mfun);
%     csc_tsd.cfg.history.cfg = cat(1,csc_tsd.cfg.history.cfg,{cfg});
% end
