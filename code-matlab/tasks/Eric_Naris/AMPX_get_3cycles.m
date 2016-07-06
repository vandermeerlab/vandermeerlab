function [all_cycles] = AMPX_get_3cycles(cfg_in, evts_in, data_in, ExpKeys)
%% Naris_get_3cycles: finds the highest amplitude peak in filtered signal
%  for each event and extracts it along with one cycle before and one cycle
%  after. Used for averaging over events for phase and CSD analyses.
%
%          Inputs:
%           - cfg [struct]: contains configuration options
%           - evts [struct]: IV struct for all the events
%           - data_tsd_remap [struct]: contains tvec, data [nChannels x samples], hdr, labels {Channels}
%           - ExpKeys [struct]: used for channel selection
%          Outputs:
%           - cycles [struct]: contains the data samples for each event.
%           -
%           -
%
% EC - 2016-05-18
%% set parameters
cfg_def.buffer = 0.2; % time (s) to include before and after gamma
cfg_def.interp_factor = 0.1; % new spacing for interpolation
cfg_def.chanlabel = diag(reshape(ExpKeys.Probe_layout,8,8)); % vertical section with no bad channels to plot for raw data.  Default for A8x8 AMPX
cfg_def.scalefactor = 600; % scale factor for raw data plot
cfg_def.Fsvid = 10; % sampling frequency for making the video
cfg_def.clims = [-800 1000]; % range of colour scale in LFP image
cfg_def.clims_csd = [-700 700]; % range of colour scale in CSD image
cfg_def.font_size = 13; % fontsize
cfg_def.t_font_size = 14; % fontsize for title
cfg_def.data_smooth = 10;
cfg_def.session_type = 'pre';
cfg_def.smooth = 10;
cfg_def.plot = 'off';
cfg_def.MPD = data_in.hdr.Fs/200;
cfg_def.peak_to_use = 'next';
cfg = ProcessConfig2(cfg_def,cfg_in);



for ievt =1:length(evts_in.tstart);
    disp(num2str(ievt))
    i1 = evts_in.tstart(ievt)*data_in.hdr.Fs;
    i2 = evts_in.tend(ievt)*data_in.hdr.Fs;
    
    % make data variable to avoid having to reload original data
    clear data2;
    data2 = data_in;
    
    evt_length = i2-i1;
    buffer = cfg.buffer*data2.hdr.Fs;
    %restrict data to section containing gamma event of interest
    if i1-buffer < 0 || i2+buffer > length(data2.tvec)
        continue
    end
    data2.tvec = data_in.tvec(i1-buffer:i2+buffer);
    for ii = 1:64
        data2.channels{ii} = data_in.channels{ExpKeys.Probe_layout(ii)}(i1-buffer:i2+buffer);
    end
    clear ii;
    
    %% filter the data BEFORE CSD!!!!!
    % FILTER THE ACTUAL DATA TO REMOVE WEIRDNESS!!!!
    data_smooth = data2;
    for ichan = 64:-1:1
        data_smooth.channels{ichan} = smooth(data2.channels{ichan}, 10);
    end
    
    data_tsd = AMPX2tsd(data_smooth); % convert to tsd
    
    %filter into band of interest
    cfg.f = cfg.freq; %filter for band of interest
    data_tsd_filt = FilterLFP(cfg, data_tsd);
    
    %% get extract the cycle
    cfg.evt_length = evt_length;
    cfg.chan = Naris_BestChan_remap(ExpKeys,'location', 'vl');
    [cycles, h_cycle] = AMPX_extract_max_cycle(cfg, data_tsd_filt); % extract the center max cycle and the one before and after it.
    peaks = {'prev', 'center', 'next'};
    
    if isempty(cycles) ==0
        for ipeak = 1:length(peaks)
            keep_idx = cycles.(peaks{ipeak}).idx; % determine which cycle to use with the cycles subfield.
            
            % keep_idx_all.(peaks{ipeak}){ievt} = keep_idx;
            
            data_tsd_filt_cycle.tvec = data_tsd_filt.tvec(keep_idx(1):keep_idx(2));
            data_tsd_filt_cycle.data = data_tsd_filt.data(:,keep_idx(1):keep_idx(2));
            %% make all the cycles the same number of points
            
            tvec_new = linspace(data_tsd_filt_cycle.tvec(1), data_tsd_filt_cycle.tvec(end), 18*4); % specific to event
            for iChan = 1:size(data_tsd_filt_cycle.data,1)
                temp(iChan, :) = interp1(data_tsd_filt_cycle.tvec, data_tsd_filt_cycle.data(iChan,:), tvec_new, 'linear');
            end
            data_tsd_filt_cycle.(peaks{ipeak}).data = temp;
            data_tsd_filt_cycle.(peaks{ipeak}).tvec = tvec_new;
            all_cycles.peaks.(peaks{ipeak}).data{ievt} = data_tsd_filt_cycle.(peaks{ipeak}).data;
            all_cycles.peaks.(peaks{ipeak}).tvec{ievt} = data_tsd_filt_cycle.(peaks{ipeak}).tvec;
        end
        
        %% append the data together to make a three cycle snapshot for the CSD
        cycle_data.data = [data_tsd_filt_cycle.prev.data(:,1:end-1), data_tsd_filt_cycle.center.data(:,1:end-1), data_tsd_filt_cycle.next.data];
        cycle_data.tvec = [data_tsd_filt_cycle.prev.tvec(:,1:end-1), data_tsd_filt_cycle.center.tvec(:,1:end-1), data_tsd_filt_cycle.next.tvec];
        %
        
        all_cycles.data{ievt} = cycle_data.data;
        all_cycles.tvec{ievt} = cycle_data.tvec;

        all_cycles.pass(ievt) = 1;
        if strcmp(cfg.plot, 'on')
            saveas(h_cycle, [datestr(date, 'yyyy-mm-dd') '/CSD/'  datestr(date, 'yyyy-mm-dd') '_' num2str(ievt) '_CSD_cycle.fig'])
        end
    else   % this is incase the event does not meet criteria.
        all_cycles.pass(ievt) = 0;
    end
    close all
end
all_cycles.cfg.chan_used = cfg.chan;
all_cycles.cfg.f = cfg.f;
all_cycles.ExpKeys = ExpKeys;

end