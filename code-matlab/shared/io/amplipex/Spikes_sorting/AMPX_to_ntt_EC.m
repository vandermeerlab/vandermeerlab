function AMPX_to_ntt_EC(cfg)
%% AMPX_to_ntt_EC: This is an expansion of the AMPX_to_ntt_JG that assumes
% that each tetrode can have spikes in both the positive and negative
% directions.  It runs through each tetrode once with positive thresholds
% and once for negative thresholds ultimately generating two ntt files
% (*TT5.ntt & *TT5_pos.ntt).
%% load data of interest
% cfg.fname = 'R060-2015-02-03-pre.dat';
fname = cfg.fname;
% cfg.cd = ['G:\Naris\' fname(1:4) '\' fname(1:15)];
cd(cfg.cd)

addpath(genpath('C:\Users\mvdm\Documents\GitHub\vandermeerlab\code-matlab\toolboxes\ums2k_02_23_2012'));

tts_to_process = cfg.tts_to_process;
cfg.ref_tt = 0;
cfg.useAverage = 0;
% cfg.thresh = [-100 -100 -100 -100]; % set this carefully
disp(['Thresholds: ' num2str(cfg.thresh)])
cfg.showData = 1;
% mkdir(['G:\Naris\' fname(1:4) '\' fname(1:end-8) '\FD'])
clear av;
for iTT = 1:length(tts_to_process)
    
    %ch = [28 26 24 22]; % tt1
    current_tt = tts_to_process(iTT);
     cfg.tts_to_process = current_tt;
    channels = tts(cfg);
    
    fprintf('Loading tetrode %d...\n',current_tt);
    data = AMPX_loadData(fname,channels.channels_array);
    
    %% OPTIONAL re-reference -- doesn't seem to make much difference when I tried this on the mPFC tt's
    if cfg.useAverage
        tic
        
        ref_channels = tts_R060(cfg.ref_tt);
        
        if ~exist('av','var')
            av = AMPX_getAverage(fname,ref_channels.channels_array); % also slow!
        end
        
        data = AMPX_reref(data,av);
        toc
    end
    
    %% convert back to AD bits for writing
    rgINT16 = (2^16)./2;
    for iC = 1:length(data.channels)
        data.channels{iC} = double(data.channels{iC}).*rgINT16; % convert to fraction of full range
        data.channels{iC} = data.channels{iC}./(data.hdr.range1_volts*10^6); % convert to microvolts
        data.channels{iC} = round(data.channels{iC}.*data.hdr.Gain); % correct for amplifier gain
    end
    
    %% resample and filter for spikes
    clear D;
    invert_waveforms = 0; % set to 1 if you have spikes bigger in the positive direction
    
    new_Fs = 32000;
    new_tvec = data.tvec(1):1./new_Fs:data.tvec(end);
    
    tic
    disp('filter')
    for ii = length(data.channels):-1:1
        % resample at 32kHz required by .ntt format
        D(ii,:) = resample(data.channels{ii},new_Fs,data.hdr.Fs);
        D(ii,:) = filter_for_spikes(D(ii,:),'Fs',new_Fs);
    end
    if invert_waveforms, D = -D; end
    toc
    
    %% take a look
    if cfg.showData
        t_end = 20; % number of seconds to plot
        
        figure(1)
        clf
        
        plot(D(:,1:t_end*new_Fs)')
    end
    
    %% detect spikes and align spikes
    for inv = 1:2
        if inv == 2; pos = 1; else pos = 0; end
        spikes = ss_default_params(new_Fs);
        spikes.params.detect_method = 'manual'; % could be 'auto'
        spikes.params.thresh = cfg.thresh;
        spikes.params.window_size = 1;
        spikes.params.cross_time = 0.4;
        if pos ==1
            spikes.params.thresh = spikes.params.thresh*-1;
        end
        tic
        disp('ss_detect')
        spikes = ss_detect({double(D')},spikes);
        toc
        
        %%
        tic
        disp('ss_align')
        spikes = ss_align(spikes);
        toc
        
        %% load test Nlx tt -- need this to have a Header to write
        load NLX_ExampleHeader;
        
        %% reformat
        nSpikes = length(spikes.spiketimes);
        
        TimestampsOUT = double(round(spikes.spiketimes*10^6));
        ScNumbersOUT = 4*ones(size(TimestampsOUT));
        CellNumbersOUT = zeros(size(TimestampsOUT));
        FeaturesOUT = zeros(8,nSpikes);
        HeaderOUT = Header;
        SamplesOUT = nan(32,4,nSpikes);
        
        for iS = 1:nSpikes
            SamplesOUT(:,:,iS) = -int16(round(spikes.waveforms(iS,:,:))*64);
        end
        
        %% write test Nlx tt
        if pos==1
            fname_out = cat(2,fname(1:end-4),'-TT',num2str(current_tt),'_pos.ntt');
        else
            fname_out = cat(2,fname(1:end-4),'-TT',num2str(current_tt),'.ntt');
        end
        Mat2NlxSpike(fname_out, 0, 1, [], [1 1 1 1 1], TimestampsOUT, ScNumbersOUT, CellNumbersOUT, FeaturesOUT, SamplesOUT, Header);
    end
end % over tetrodes
end
