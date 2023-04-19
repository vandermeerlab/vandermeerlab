%% Script to calculate timeresolved PPC for cells that have significant phase locking

% setup parameters
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items

% Load the significantly phase_locked cells
load('D:\RandomVstrAnalysis\final_results\significantly_phase_locked.mat');

% Cycle through FSIs first
for iC = 1:size(sig3_fsi,1)
    this_sig = sig3_fsi{iC,2};
    % Skip this cell if it is not significantly phase locked to any
    % frequency band
    if ~any(this_sig)
        continue
    end
    label = sig3_fsi{iC,1};
    % generate correct path on the basis of the label
    toks = strsplit(label, '-');
    path = strcat('E:\ADRLabData\', toks{1}, '\', strjoin(toks(1:4),'-'));
    cd(path);
    
    % Load CSC, ExpKeys and the cell
    LoadExpKeys;
    
    if isfield(ExpKeys,'goodGamma_vStr')
        cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
    elseif isfield(ExpKeys, 'goodGamma')
        cfg = []; cfg.fc = ExpKeys.goodGamma;
    else
        error('Couldn''t find LFP field name.');
    end
    
    csc = LoadCSC(cfg); csc.data = csc.data-nanmean(csc.data); % could locdetrend to improve STA estimate
    % Note that Fieldtrip will interpoalate data to avoid gaps. Include
    % sanity tests to ensure STA/STS segments don't include these gaps.
    ft_csc = ft_read_neuralynx_interp(cfg.fc);
    
    % hacky code to separate out only ontrack-data
    temp_tvec = [0:length(ft_csc.time{1})-1];
    temp_offset = (double(ft_csc.hdr.LastTimeStamp)/1e6 - double(ft_csc.hdr.FirstTimeStamp)/1e6)/(length(temp_tvec) - 1);
    temp_tvec = temp_tvec * temp_offset;
     
    % Modify ft_csc
    ft_csc.time{1} = temp_tvec;
    ft_csc.fsample = 1/temp_offset;
    
    % Ensure that the new data doesn't have any Nans in it
    temp_tvec = temp_tvec + double(ft_csc.hdr.FirstTimeStamp)/1e6;
    temp_start = nearest_idx3(ExpKeys.TimeOnTrack, temp_tvec);
    temp_end = nearest_idx3(ExpKeys.TimeOffTrack, temp_tvec);
    cfg_onTrack.begsample = temp_start;
    cfg_onTrack.endsample = temp_end;
    data_onTrack = ft_redefinetrial(cfg_onTrack, ft_csc);
    
    if ~isempty(find(isnan(data_onTrack.time{1}),1))
        warning('On Track data for %s has gaps', cfg.fc{1});
    end

    cfg = []; cfg.fc = {strcat(label, '.t')};
    sd.S = LoadSpikes(cfg);
    % Restrict spikes to only OnTrack
    sd.S = restrict(sd.S, iv(ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack));
    % Read the spike files into field trip format
    sd.S.ft_spikes = ft_read_spike(sd.S.label{1});
    
    % Calculate and save STA
    cfg_ft.timwin = [-0.5 0.5];
    cfg_ft.spikechannel = sd.S.ft_spikes(1).label{1};
    cfg_ft.channel = ft_csc.label(1);
    this_data = ft_appendspike([], ft_csc, sd.S.ft_spikes(1));

    % restrict spikes to a timeWindow of +/-5 seconds around the first reward fire time  
    rt1 = getRewardTimes();
    rt1 = sort(rt1);
    rt1 = rt1(rt1 > ExpKeys.TimeOnTrack);
    % Sometimes (in R117-2007-06-12, for instance) getRewardTimes() returns
    % times that are spaced out less than 5 sec apart (possibly erroneus). 
    % Getting rid of such reward times to maintain consistency
    rt_dif = diff(rt1);
    rt_dif = find(rt_dif <= 5);
    valid_rt1 = true(length(rt1),1);
    for i = 1:length(rt_dif)
        valid_rt1(rt_dif(i)) = false;
        valid_rt1(rt_dif(i)+1) = false;
    end
    rt1 = rt1(valid_rt1);
   
    
    % For near reward_trials  
    w_start = rt1 - 5;
    w_end =  rt1 + 5;
    % If the last trial time exceeds end of experiment, discard that trial
    if w_end(end) < ExpKeys.TimeOffTrack
        w_end(end) = []; w_start(end) = [];
    end
    rt_iv = iv(w_start, w_end);

    % Break down data into near trials
    temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;     
    temp_start = nearest_idx3(rt_iv.tstart, temp_tvec);
    temp_end = nearest_idx3(rt_iv.tend, temp_tvec);
    cfg_near_trials.trl = [temp_start, temp_end, zeros(size(temp_start))];
    near_data = ft_redefinetrial(cfg_near_trials, this_data);
    this_ft_spk_count = 0;
    near_trialwise_spk_count = zeros(1, length(near_data.trial));
    for iT = 1:length(near_data.trial)
        near_trialwise_spk_count(iT) = sum(near_data.trial{iT}(2,:));
        this_ft_spk_count = this_ft_spk_count + near_trialwise_spk_count(iT);
    end

    % Sanity check to ensure that redefine trial has the same number of
    % spikes as restrict()
    this_spk_count = length(restrict(sd.S, rt_iv).t{1});

    % Extract All Spike IDs
    trial_wise_spike = cell(1, length(near_data.trial));
    last_spk_ct = 0;
    for iT = 1:length(near_data.trial)
        this_spk_ct = near_trialwise_spk_count(iT);
        trial_wise_spike{iT} = last_spk_ct + 1 : last_spk_ct + this_spk_ct;
        last_spk_ct = last_spk_ct + this_spk_ct;
    end
    
    % Calculate STS for all spikes
    cfg_sts.method = 'mtmconvol';
    cfg_sts.foi = 1:1:100;
    cfg_sts.t_ftimwin = 5./cfg_sts.foi; % 5 cycles per frequency
    cfg_sts.taper = 'hanning';
    cfg_sts.spikechannel =  sd.S.ft_spikes(1).label{1};
    cfg_sts.channel = near_data.label{1};
    this_sts = ft_spiketriggeredspectrum(cfg_sts, near_data);
    this_flag = false;
    % Display warning to show that there were Nans in this calculation
    if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
        this_flag = true;
        warning('Cell %s has nans in its STS',sd.S.label{iC});
    end
    this_freq = this_sts.freq;
    this_sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
    flag_nansts = this_flag;
    
    % Calculate and save trialwise_ppc and use subsampled measures
    % for near Reward Trials. However, if an FSI has fewer spikes
    % than at least one co-recorded MSN, do NOT
    % subsample and include all the spikes
    
    % Calculate non time-resolved PPC
    cfg_ppc               = [];
    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
    cfg_ppc.spikechannel  = this_sts.label;
    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
    cfg_ppc.avgoverchan   = 'weighted';
    this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc, this_sts);

    %% Calculate time-resolved PPC
    cfg_ppc2 = cfg_ppc;
    cfg_ppc2.latency = 'maxperiod';
    cfg_ppc2.timwin = 2;
    cfg_ppc2.winstepsize = 0.1;
    this_tr_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc2, this_sts);
    mean_tr_fr = squeeze(this_tr_ppc.nspikes)/size(near_data.trial,2);
    this_bands = find(this_sig);
    tp = length(this_bands) + 1;
    fig = figure('WindowState', 'maximized');
    ax = subplot(tp, 1, 1);
    plot(this_tr_ppc.time, mean_tr_fr, 'black')
    ax.XAxis.Label.String = 'Time';
    ax.XAxis.TickLabels = {'-5', '-3', '-1', '1', '3', '5'};
    ax.YAxis.Label.String = 'Inst. F.R';
    ax.XAxis.Label.FontSize = 12;
    ax.YAxis.Label.FontSize = 12;
    ax.Title.Interpreter = 'none';
    ax.Title.String = label;
    for iB = 1:length(this_bands)
        f_idx = find(round(this_freq) >= f_list{this_bands(iB)}(1) & ...
            round(this_freq) <= f_list{this_bands(iB)}(2));
        f_ppc = squeeze(this_tr_ppc.ppc0);
        f_ppc = mean(f_ppc(f_idx,:));
        ax = subplot(tp,1,iB+1);
        plot(this_tr_ppc.time, f_ppc, 'Color', c_list{this_bands(iB)})
        ax.XAxis.Label.String = 'Time';
        ax.XAxis.TickLabels = {'-5', '-3', '-1', '1', '3', '5'};
        ax.YAxis.Label.String = 'PPC';
        ax.YAxis.Exponent = 0;
        ax.Title.String = sprintf("FSI %d Hz - %d Hz", f_list{this_bands(iB)}(1), ...
            f_list{this_bands(iB)}(2));
        ax.Title.Color = c_list{this_bands(iB)};
        ax.XAxis.Label.FontSize = 12;
        ax.YAxis.Label.FontSize = 12;
    end
    WriteFig(fig, strcat('D:\RandomVstrAnalysis\PPCbyFR\FSI_',label,'_FRbyPPC'),1)
    close;
end

% Cycle through MSNs next
for iC = 1:size(sig3_msn,1)

end




