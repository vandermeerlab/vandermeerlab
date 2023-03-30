%% Script to generate PPC as a function of binned firing-rate (now split into bins such that spikes are more or less equally distributed)

% setup
clear;
%load significant cells
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
load('./significance.mat');
num_subs = 1000;
nbins = length(z_edges-1); % Number of bins
nshufs = 1000; % number of shuffles to calculate STS/PPC Significance
spk_dt = 0.025; % interspike interval for surrogate spike train used for spike-triggered spectrum pool
visited  = dictionary('dummy',{-1}); % Use this to avoid creating fake spike-triggered spectra again

% Cycle through MSNs first
for iC = 1:length(sig_msn)
    label = sig_msn{iC,1};
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
    
    % Block of code to divide recordings session into trials
    
    % restrict spikes to a timeWindow of +/-5 seconds around the reward  
    rt1 = getRewardTimes();
    rt1 = rt1(rt1 > ExpKeys.TimeOnTrack);
    rt2 = getRewardTimes2();
    rt2 = rt2(rt2 > ExpKeys.TimeOnTrack);
    % Sometimes (in R117-2007-06-12, for instance) getRewardTimes() returns
    % times that are spaced out less than 5 sec apart (possibly erroneus). 
    % Getting rid of such reward times to maintain consistency
    rt_dif = diff(rt1);
    rt_dif = find(rt_dif <= 5);
    valid_rt1 = true(length(rt1),1);
    valid_rt2 = true(length(rt2),1);
    for i = 1:length(rt_dif)
    valid_rt1(rt_dif(i)) = false;
    valid_rt1(rt_dif(i)+1) = false;
    valid_rt2(rt2 >= rt1(rt_dif(i)) & rt2 <= rt1(rt_dif(i)+2)) = false;
    end
    % Sometimes (in R119-2007-07-05, for instance) getRewardTimes2() returns
    % times that are spaced out less than 5 sec apart (possibly erroneus). 
    % Getting rid of such reward times to maintain consistency
    rt_dif = diff(rt2);
    rt_dif = find(rt_dif <= 5);
    for i = 1:length(rt_dif)
    valid_rt2(rt_dif(i)) = false;
    valid_rt2(rt_dif(i)+1) = false;
    valid_rt1(rt1 >= rt2(rt_dif(i)-1) & rt1 <= rt2(rt_dif(i)+1)) = false;
    end
    rt1 = rt1(valid_rt1);
    rt2 = rt2(valid_rt2);
    % Sometimes (in R117-2007-06-17, for instance) the second reward is
    % triggered but not the first one in the last trial
    if length(rt1) ~= length(rt2)
    rt1 = rt1(1:end-1);
    end
    % Sanity check to make sure that rt2 is always triggered after rt1
    keep = (rt1 <= rt2);
    rt1 = rt1(keep);
    rt2 = rt2(keep);
    
    % For near reward_trials  
    w_start = rt1 - 5;
    w_end =  rt2 + 5;
    % Last trial time shouldn't exceed Experiment end time
    w_end(end) = min(w_end(end), ExpKeys.TimeOffTrack);
    % Sorting makes it wonky in some cases (in R119-2007-07-06),so 
    % only keep trials that are not outliers
    keep = ~isoutlier(w_end - w_start, 'median');
    w_start = w_start(keep);
    w_end = w_end(keep);
    rt_iv = iv(w_start, w_end);
    
    % Break down data into near trials
    temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;     
    temp_start = nearest_idx3(rt_iv.tstart, temp_tvec);
    temp_end = nearest_idx3(rt_iv.tend, temp_tvec);
    cfg_near_trials.trl = [temp_start, temp_end, zeros(size(temp_start))];
    near_data = ft_redefinetrial(cfg_near_trials, this_data);
    % Sanity check to ensure that redefine trial has the same number of
    % spikes as restrict()
    this_flag = false;
    spk_count1 = length(restrict(sd.S, rt_iv).t{1});
    spk_count2 = 0;
    for iT = 1:length(near_data.trial)
    spk_count2 = spk_count2 + sum(near_data.trial{iT}(2,:)); 
    end
    od.spk_count = spk_count2;
    od.flag_tooFewSpikes = false;
    if spk_count1 ~= spk_count2
        this_flag = true;
        warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in near-Reward trials', ...
        spk_count1, spk_count2)
    end

    % Calculate and save STA
    this_sta = ft_spiketriggeredaverage(cfg_ft, near_data);
    od.sta_time = this_sta.time;
    od.unsampled_sta = this_sta.avg(:,:)';
    od.flag_unequalSpikes = this_flag;
    
    % Calculate and save STS
    cfg_sts.method = 'mtmconvol';
    cfg_sts.foi = 1:1:100;
    cfg_sts.t_ftimwin = 5./cfg_sts.foi;
    cfg_sts.taper = 'hanning';
    cfg_sts.spikechannel =  sd.S.ft_spikes(1).label{1};
    cfg_sts.channel = near_data.label{1};
    this_sts = ft_spiketriggeredspectrum(cfg_sts, near_data);
    this_flag = false;
    % Display warning to show that there were Nans in this calculation
    if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
        this_flag = true;
        warning('Cell %s has nans in its STS',sd.S.label{1});
    end
    od.freqs = this_sts.freq;
    od.unsampled_sts = nanmean(sq(abs(this_sts.fourierspctrm{1})));
    od.flag_nansts = this_flag;

    % Configure ppc
    cfg_ppc               = [];
    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
    cfg_ppc.spikechannel  = this_sts.label;
    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
    cfg_ppc.avgoverchan   = 'weighted';
    cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window

    if(~isKey(visited, path))
        % If visiting this path for the first time, create session wide STS-pool
        cfg_f.begsample = cfg_onTrack.begsample;
        cfg_f.endsample = cfg_onTrack.endsample;
        f_spike = sd.S.ft_spikes;
        f_spike.timestamp{1} = double(ft_csc.hdr.FirstTimeStamp) +  ...
            10^6*(ft_csc.time{1}(1):spk_dt:ft_csc.time{1}(end));
        f_data = ft_appendspike([], ft_csc,f_spike);
        f_data = ft_redefinetrial(cfg_f,f_data);
        fprintf("Number of fake spikes for session %s is %d\n", path, sum(f_data.trial{1}(2,:)));
        cfg_f = [];
        cfg_f.method = 'mtmconvol';
        cfg_f.foi = 1:1:100;
        cfg_f.t_ftimwin = 5./cfg_sts.foi;
        cfg_f.taper = 'hanning';
        cfg_f.spikechannel =  f_spike.label{1};
        cfg_sts.channel = f_data.label{1};
        f_sts = ft_spiketriggeredspectrum(cfg_f, f_data);
        visited(path) = {f_sts};
        clear f_data f_sts cfg_f f_spike
    else
        fprintf('This saved you time!\n')
    end

    % do shuffles for significance_testing
    od.nshufs = nshufs;
    pool_sts = visited(path);
    pool_sts = pool_sts{1};
    shuf_sts = zeros(nshufs, length(od.unsampled_sts));
    shuf_ppc = zeros(nshufs, length(od.unsampled_sts));
    for iShuf = 1:nshufs
        pool_count = length(pool_sts.time{1});
        keep = randperm(pool_count);
        keep = keep(1:od.spk_count);
        this_shuf = pool_sts;
        this_shuf.label{1} = label; % Because the pooled STS might have a differnt spike channel label
        this_shuf.fourierspctrm{1} = pool_sts.fourierspctrm{1}(keep,:,:);
        this_shuf.time{1} = pool_sts.time{1}(keep,:);
        this_shuf.trial{1} = pool_sts.trial{1}(keep,:);
        shuf_sts(iShuf,:) = nanmean(sq(abs(this_shuf.fourierspctrm{1})));
        this_shuf_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, this_shuf);
        shuf_ppc(iShuf,:) = this_shuf_ppc.ppc0;
    end
    od.shuf_sts = shuf_sts;
    od.shuf_ppc = shuf_ppc;

    % Get mfr to bin trials
    tcount = length(rt_iv.tstart);
    mfr = zeros(tcount, 1);
    all_tspikes = cell(1,tcount);
    for iT = 1:tcount
        all_tspikes{iT} = sum(near_data.trial{iT}(2,:));
        mfr(iT) = all_tspikes{iT}/near_data.time{iT}(end); 
    end
    od.mfr = mfr;

    % Bin trials by z_scored firing rate and then calculate PPC for each bin
    
    % get rid of trials with 0 spikes but then put them back together
    % later with the 0 spike trials having an absurd bin value (like -1)
    nz_trials = mfr > 0;
    fr_bin = zeros(size(mfr));
    fr_bin(~nz_trials) = -1;
    nz_mfr = mfr(nz_trials);
    tw_fr = zscore(nz_mfr);
    [bcount, edges, nz_bin] = histcounts(tw_fr, z_edges);
    % Redistributing bin values to take care of zero-spike-trials
    inZ = 1;
    for iZ = 1:length(fr_bin)
        if nz_trials(iZ)
            fr_bin(iZ) = nz_bin(inZ);
            inZ = inZ + 1;
        end
    end
    [od.z_bcount, od.z_edges, od.z_fr_bin] = deal(bcount, edges, fr_bin);
    od.z_binned_mean = mean(nz_mfr);
    od.z_binned_sd = std(nz_mfr);
    od.z_binned_spk_count = zeros(length(bcount),1);
    od.z_binned_sts = nan(length(bcount), length(this_sts.freq));
    od.z_binned_ppc = nan(length(bcount), length(this_sts.freq));
    trl_idx = this_sts.trial{1};        
    bin_idx = zeros(size(trl_idx));

    for iB = 1:length(bcount)
        % get corresponding fr_bins for spikes
        this_trials = find(fr_bin == iB);
        for iT = 1:length(this_trials)
            bin_idx(trl_idx == this_trials(iT)) = iB;
        end
        bin_sts = this_sts;
        bin_sts.fourierspctrm{1} = bin_sts.fourierspctrm{1}(bin_idx == iB,:,:);
        bin_sts.time{1} = bin_sts.time{1}(bin_idx == iB,:);
        bin_sts.trial{1} = bin_sts.trial{1}(bin_idx == iB,:);
        od.z_binned_sts(iB,:) = nanmean(sq(abs(bin_sts.fourierspctrm{1})));
        bin_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, bin_sts);
        od.z_binned_ppc(iB,:) = bin_ppc.ppc0;
        od.z_binned_spk_count(iB) = sum(bin_idx == iB);
    end

    % Bin trials by normalizing firing rate between 0-1 and then calculate PPC for each bin
 
    % get rid of trials with 0 spikes but then put them back together
    % later with the 0 spike trials having an absurd bin value (like -1)
    nz_trials = mfr > 0;
    fr_bin = zeros(size(mfr));
    fr_bin(~nz_trials) = -1;
    nz_mfr = mfr(nz_trials);
    tw_fr = (nz_mfr-min(nz_mfr))/(max(nz_mfr) - min(nz_mfr));
    [bcount, edges, nz_bin] = histcounts(tw_fr, n_edges);
    % Redistributing bin values to take care of zero-spike-trials
    inZ = 1;
    for iZ = 1:length(fr_bin)
        if nz_trials(iZ)
            fr_bin(iZ) = nz_bin(inZ);
            inZ = inZ + 1;
        end
    end
    [od.n_bcount, od.n_edges, od.n_fr_bin] = deal(bcount, edges, fr_bin);
    od.n_binned_spk_count = zeros(length(bcount),1);
    od.n_binned_sts = nan(length(bcount), length(this_sts.freq));
    od.n_binned_ppc = nan(length(bcount), length(this_sts.freq));
    trl_idx = this_sts.trial{1};        
    bin_idx = zeros(size(trl_idx));

    for iB = 1:length(bcount)
        % get corresponding fr_bins for spikes
        this_trials = find(fr_bin == iB);
        for iT = 1:length(this_trials)
            bin_idx(trl_idx == this_trials(iT)) = iB;
        end
        bin_sts = this_sts;
        bin_sts.fourierspctrm{1} = bin_sts.fourierspctrm{1}(bin_idx == iB,:,:);
        bin_sts.time{1} = bin_sts.time{1}(bin_idx == iB,:);
        bin_sts.trial{1} = bin_sts.trial{1}(bin_idx == iB,:);
        od.n_binned_sts(iB,:) = nanmean(sq(abs(bin_sts.fourierspctrm{1})));
        bin_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, bin_sts);
        od.n_binned_ppc(iB,:) = bin_ppc.ppc0;
        od.n_binned_spk_count(iB) = sum(bin_idx == iB);
    end

    % Save values
    save(strcat('D:\RandomVstrAnalysis\temp\',label, '_MSN_ppcbyfr.mat'), 'od')
end

% Cycle through FSIs next
for iC = 1:length(sig_fsi)  
    label = sig_fsi{iC,1};
    toks = strsplit(label, '-');

    % Load the results file to get the msn spike_count distribution
    full_res = load(strcat('D:\RandomVstrAnalysis\final_results\', strjoin(toks(1:4),'-'),'_ft_spec.mat'));
    msn_dist = full_res.od.msn_near_dist;

    % generate correct path on the basis of the label
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
    
    % Block of code to divide recordings session into trials
    
    % restrict spikes to a timeWindow of +/-5 seconds around the reward  
    rt1 = getRewardTimes();
    rt1 = rt1(rt1 > ExpKeys.TimeOnTrack);
    rt2 = getRewardTimes2();
    rt2 = rt2(rt2 > ExpKeys.TimeOnTrack);
    % Sometimes (in R117-2007-06-12, for instance) getRewardTimes() returns
    % times that are spaced out less than 5 sec apart (possibly erroneus). 
    % Getting rid of such reward times to maintain consistency
    rt_dif = diff(rt1);
    rt_dif = find(rt_dif <= 5);
    valid_rt1 = true(length(rt1),1);
    valid_rt2 = true(length(rt2),1);
    for i = 1:length(rt_dif)
    valid_rt1(rt_dif(i)) = false;
    valid_rt1(rt_dif(i)+1) = false;
    valid_rt2(rt2 >= rt1(rt_dif(i)) & rt2 <= rt1(rt_dif(i)+2)) = false;
    end
    % Sometimes (in R119-2007-07-05, for instance) getRewardTimes2() returns
    % times that are spaced out less than 5 sec apart (possibly erroneus). 
    % Getting rid of such reward times to maintain consistency
    rt_dif = diff(rt2);
    rt_dif = find(rt_dif <= 5);
    for i = 1:length(rt_dif)
    valid_rt2(rt_dif(i)) = false;
    valid_rt2(rt_dif(i)+1) = false;
    valid_rt1(rt1 >= rt2(rt_dif(i)-1) & rt1 <= rt2(rt_dif(i)+1)) = false;
    end
    rt1 = rt1(valid_rt1);
    rt2 = rt2(valid_rt2);
    % Sometimes (in R117-2007-06-17, for instance) the second reward is
    % triggered but not the first one in the last trial
    if length(rt1) ~= length(rt2)
    rt1 = rt1(1:end-1);
    end
    % Sanity check to make sure that rt2 is always triggered after rt1
    keep = (rt1 <= rt2);
    rt1 = rt1(keep);
    rt2 = rt2(keep);
    
    % For near reward_trials  
    w_start = rt1 - 5;
    w_end =  rt2 + 5;
    % Last trial time shouldn't exceed Experiment end time
    w_end(end) = min(w_end(end), ExpKeys.TimeOffTrack);
    % Sorting makes it wonky in some cases (in R119-2007-07-06),so 
    % only keep trials that are not outliers
    keep = ~isoutlier(w_end - w_start, 'median');
    w_start = w_start(keep);
    w_end = w_end(keep);
    rt_iv = iv(w_start, w_end);
    
    % Break down data into near trials
    temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;     
    temp_start = nearest_idx3(rt_iv.tstart, temp_tvec);
    temp_end = nearest_idx3(rt_iv.tend, temp_tvec);
    cfg_near_trials.trl = [temp_start, temp_end, zeros(size(temp_start))];
    near_data = ft_redefinetrial(cfg_near_trials, this_data);
    % Sanity check to ensure that redefine trial has the same number of
    % spikes as restrict()
    this_flag = false;
    spk_count1 = length(restrict(sd.S, rt_iv).t{1});
    spk_count2 = 0;
    for iT = 1:length(near_data.trial)
    spk_count2 = spk_count2 + sum(near_data.trial{iT}(2,:)); 
    end
    od.spk_count = spk_count2;
    od.flag_tooFewSpikes = false;
    if spk_count1 ~= spk_count2
        this_flag = true;
        warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in near-Reward trials', ...
        spk_count1, spk_count2)
    end

    % Calculate and save STA
    this_sta = ft_spiketriggeredaverage(cfg_ft, near_data);
    od.sta_time = this_sta.time;
    od.unsampled_sta = this_sta.avg(:,:)';
    od.flag_unequalSpikes = this_flag;
    
    % Calculate and save STS
    cfg_sts.method = 'mtmconvol';
    cfg_sts.foi = 1:1:100;
    cfg_sts.t_ftimwin = 5./cfg_sts.foi;
    cfg_sts.taper = 'hanning';
    cfg_sts.spikechannel =  sd.S.ft_spikes(1).label{1};
    cfg_sts.channel = near_data.label{1};
    this_sts = ft_spiketriggeredspectrum(cfg_sts, near_data);
    this_flag = false;
    % Display warning to show that there were Nans in this calculation
    if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
        this_flag = true;
        warning('Cell %s has nans in its STS',sd.S.label{1});
    end
    od.freqs = this_sts.freq;
    od.unsampled_sts = nanmean(sq(abs(this_sts.fourierspctrm{1})));
    od.flag_nansts = this_flag;

    % Configure ppc
    cfg_ppc               = [];
    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
    cfg_ppc.spikechannel  = this_sts.label;
    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
    cfg_ppc.avgoverchan   = 'weighted';
    cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window

    if(~isKey(visited, path))
        % If visiting this path for the first time, create session wide STS-pool
        cfg_f.begsample = cfg_onTrack.begsample;
        cfg_f.endsample = cfg_onTrack.endsample;
        f_spike = sd.S.ft_spikes;
        f_spike.timestamp{1} = double(ft_csc.hdr.FirstTimeStamp) +  ...
            10^6*(ft_csc.time{1}(1):spk_dt:ft_csc.time{1}(end));
        f_data = ft_appendspike([], ft_csc,f_spike);
        f_data = ft_redefinetrial(cfg_f,f_data);
        fprintf("Number of fake spikes for session %s is %d\n", path, sum(f_data.trial{1}(2,:)));
        cfg_f = [];
        cfg_f.method = 'mtmconvol';
        cfg_f.foi = 1:1:100;
        cfg_f.t_ftimwin = 5./cfg_sts.foi;
        cfg_f.taper = 'hanning';
        cfg_f.spikechannel =  f_spike.label{1};
        cfg_sts.channel = f_data.label{1};
        f_sts = ft_spiketriggeredspectrum(cfg_f, f_data);
        visited(path) = {f_sts};
        clear f_data f_sts cfg_f f_spike
    else
        fprintf('This saved you time!')
    end

    % do shuffles for significance_testing
    od.nshufs = nshufs;
    pool_sts = visited(path);
    pool_sts = pool_sts{1};
    shuf_sts = zeros(nshufs, length(od.unsampled_sts));
    shuf_ppc = zeros(nshufs, length(od.unsampled_sts));
    for iShuf = 1:nshufs
        pool_count = length(pool_sts.time{1});
        keep = randperm(pool_count);
        keep = keep(1:od.spk_count);
        this_shuf = pool_sts;
        this_shuf.label{1} = label; % Because the pooled STS might have a differnt spike channel label
        this_shuf.fourierspctrm{1} = pool_sts.fourierspctrm{1}(keep,:,:);
        this_shuf.time{1} = pool_sts.time{1}(keep,:);
        this_shuf.trial{1} = pool_sts.trial{1}(keep,:);
        shuf_sts(iShuf,:) = nanmean(sq(abs(this_shuf.fourierspctrm{1})));
        this_shuf_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, this_shuf);
        shuf_ppc(iShuf,:) = this_shuf_ppc.ppc0;
    end
    od.shuf_sts = shuf_sts;
    od.shuf_ppc = shuf_ppc;

    % Get mfr to bin trials
    tcount = length(rt_iv.tstart);
    mfr = zeros(tcount, 1);
    all_tspikes = cell(1,tcount);
    for iT = 1:tcount
        all_tspikes{iT} = sum(near_data.trial{iT}(2,:));
        mfr(iT) = all_tspikes{iT}/near_data.time{iT}(end); 
    end
    od.mfr = mfr;

    % Check if subsampling is needed or not
    all_spikes = sum(cell2mat(all_tspikes));
    if all_spikes < max(msn_dist)
        to_sub = false;
    else
        to_sub = true;
    end
    
    % Bin trials by z_scored firing rate and then calculate PPC for each bin

    % get rid of trials with 0 spikes but then put them back together
    % later with the 0 spike trials having an absurd bin value (like -1)
    nz_trials = mfr > 0;
    fr_bin = zeros(size(mfr));
    fr_bin(~nz_trials) = -1;
    nz_mfr = mfr(nz_trials);
    tw_fr = zscore(nz_mfr);
    [bcount, edges, nz_bin] = histcounts(tw_fr, z_edges);
    % Redistributing bin values to take care of zero-spike-trials
    inZ = 1;
    for iZ = 1:length(fr_bin)
        if nz_trials(iZ)
            fr_bin(iZ) = nz_bin(inZ);
            inZ = inZ + 1;
        end
    end
    [od.z_bcount, od.z_edges, od.z_fr_bin] = deal(bcount, edges, fr_bin);
    od.z_binned_mean = mean(nz_mfr);
    od.z_binned_sd = std(nz_mfr);
    od.z_binned_spk_count = zeros(length(bcount),1);
    od.z_binned_sts = nan(length(bcount), length(this_sts.freq));
    od.z_binned_ppc = nan(length(bcount), length(this_sts.freq));
    trl_idx = this_sts.trial{1};        
    bin_idx = zeros(size(trl_idx));
    msn_dist = round(msn_dist/length(bcount)); % Converting MSN Dist to an expected per bin count
      
    for iB = 1:length(bcount)
        % get corresponding fr_bins for spikes
        this_trials = find(fr_bin == iB);
        for iT = 1:length(this_trials)
            bin_idx(trl_idx == this_trials(iT)) = iB;
        end
        bin_spks = find(bin_idx == iB);
        if  ~isempty(bin_spks)
            if to_sub && (length(bin_spks) > round(max(msn_dist)))
                temp_ppc = zeros(num_subs,length(this_sts.freq));
                temp_sts = zeros(num_subs,length(this_sts.freq));
                for iSub = 1:num_subs
                    sub_count = msn_dist(randi(length(msn_dist)));
                    sub_idx = sort(randperm(length(bin_spks), sub_count));
                    sub_spks = bin_spks(sub_idx);
                    bin_sts = this_sts;
                    bin_sts.fourierspctrm{1} = bin_sts.fourierspctrm{1}(sub_spks,:,:);
                    bin_sts.time{1} = bin_sts.time{1}(sub_spks,:);
                    bin_sts.trial{1} = bin_sts.trial{1}(sub_spks,:);
                    temp_sts(iSub,:) = nanmean(sq(abs(bin_sts.fourierspctrm{1})));
                    bin_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, bin_sts);
                    temp_ppc(iSub,:) = bin_ppc.ppc0;
                end
                    od.z_binned_sts(iB,:) = mean(temp_sts);
                    od.z_binned_ppc(iB,:) = mean(temp_ppc);
                    od.z_binned_spk_count(iB) = length(bin_spks);
            else
                bin_sts = this_sts;
                bin_sts.fourierspctrm{1} = bin_sts.fourierspctrm{1}(bin_spks,:,:);
                bin_sts.time{1} = bin_sts.time{1}(bin_spks,:);
                bin_sts.trial{1} = bin_sts.trial{1}(bin_spks,:);
                od.z_binned_sts(iB,:) = nanmean(sq(abs(bin_sts.fourierspctrm{1})));
                bin_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, bin_sts);
                od.z_binned_ppc(iB,:) = bin_ppc.ppc0;
                od.z_binned_spk_count(iB) = length(bin_spks);
            end
        end
    end

    % Bin trials by normalizing firing rate between 0-1 and then calculate PPC for each bin

    % get rid of trials with 0 spikes but then put them back together
    % later with the 0 spike trials having an absurd bin value (like -1)
    nz_trials = mfr > 0;
    fr_bin = zeros(size(mfr));
    fr_bin(~nz_trials) = -1;
    nz_mfr = mfr(nz_trials);
    tw_fr = (nz_mfr - min(nz_mfr))/(max(nz_mfr) - min(nz_mfr));
    [bcount, edges, nz_bin] = histcounts(tw_fr, n_edges);
    % Redistributing bin values to take care of zero-spike-trials
    inZ = 1;
    for iZ = 1:length(fr_bin)
        if nz_trials(iZ)
            fr_bin(iZ) = nz_bin(inZ);
            inZ = inZ + 1;
        end
    end
    [od.n_bcount, od.n_edges, od.n_fr_bin] = deal(bcount, edges, fr_bin);
    od.n_binned_spk_count = zeros(length(bcount),1);
    od.n_binned_sts = nan(length(bcount), length(this_sts.freq));
    od.n_binned_ppc = nan(length(bcount), length(this_sts.freq));
    trl_idx = this_sts.trial{1};        
    bin_idx = zeros(size(trl_idx));
    msn_dist = round(msn_dist/length(bcount)); % Converting MSN Dist to an expected per bin count
      
    for iB = 1:length(bcount)
        % get corresponding fr_bins for spikes
        this_trials = find(fr_bin == iB);
        for iT = 1:length(this_trials)
            bin_idx(trl_idx == this_trials(iT)) = iB;
        end
        bin_spks = find(bin_idx == iB);
        if  ~isempty(bin_spks)
            if to_sub && (length(bin_spks) > round(max(msn_dist)))
                temp_ppc = zeros(num_subs,length(this_sts.freq));
                temp_sts = zeros(num_subs,length(this_sts.freq));
                for iSub = 1:num_subs
                    sub_count = msn_dist(randi(length(msn_dist)));
                    sub_idx = sort(randperm(length(bin_spks), sub_count));
                    sub_spks = bin_spks(sub_idx);
                    bin_sts = this_sts;
                    bin_sts.fourierspctrm{1} = bin_sts.fourierspctrm{1}(sub_spks,:,:);
                    bin_sts.time{1} = bin_sts.time{1}(sub_spks,:);
                    bin_sts.trial{1} = bin_sts.trial{1}(sub_spks,:);
                    temp_sts(iSub,:) = nanmean(sq(abs(bin_sts.fourierspctrm{1})));
                    bin_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, bin_sts);
                    temp_ppc(iSub,:) = bin_ppc.ppc0;
                end
                    od.n_binned_sts(iB,:) = mean(temp_sts);
                    od.n_binned_ppc(iB,:) = mean(temp_ppc);
                    od.n_binned_spk_count(iB) = length(bin_spks);
            else
                bin_sts = this_sts;
                bin_sts.fourierspctrm{1} = bin_sts.fourierspctrm{1}(bin_spks,:,:);
                bin_sts.time{1} = bin_sts.time{1}(bin_spks,:);
                bin_sts.trial{1} = bin_sts.trial{1}(bin_spks,:);
                od.n_binned_sts(iB,:) = nanmean(sq(abs(bin_sts.fourierspctrm{1})));
                bin_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, bin_sts);
                od.n_binned_ppc(iB,:) = bin_ppc.ppc0;
                od.n_binned_spk_count(iB) = length(bin_spks);
            end
        end
    end
    % Save values
    save(strcat('D:\RandomVstrAnalysis\temp\',label, '_FSI_ppcbyfr.mat'), 'od')
end

%% Test plotting

clear; 
cd('D:\RandomVstrAnalysis\PPCbyFR\');
paths = FindFiles('*FSI_ppcbyfr.mat');
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
min_spikes = 50;
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
edges = [-100,-2,-1,0,1,2,100];
[max_out, mean_out] = deal(cell(length(paths),length(f_list)));
for iC = 1:length(paths)
    load(paths{iC});
    for iF = 1:length(f_list)
        f_idx = find(round(od.freqs) >= f_list{iF}(1) & round(od.freqs) <= f_list{iF}(2));
        [binned_max, binned_mean] = deal(zeros(1,size(od.binned_ppc,1)));
        for iB = 1:size(od.binned_ppc,1)
            if od.binned_spk_count(iB) < min_spikes
                binned_max(iB) = NaN;
                binned_mean(iB) = NaN;
            else
                binned_max(iB) = max(od.binned_ppc(iB,f_idx));
                binned_mean(iB) = mean(od.binned_ppc(iB,f_idx));
            end
        end
        % Convert Inf to Nan
        binned_mean(~isfinite(binned_mean)) = NaN;
        binned_max(~isfinite(binned_max)) = NaN;
        % normalize it
        isnum = ~isnan(binned_mean);
        binned_mean(isnum) = (binned_mean(isnum) - min(binned_mean(isnum))) /...
            (max(binned_mean(isnum) - min(binned_mean(isnum))));
        isnum = ~isnan(binned_max);
        binned_max(isnum) = (binned_max(isnum) - min(binned_max(isnum))) /...
            (max(binned_max(isnum) - min(binned_max(isnum))));
        max_out{iC,iF} = binned_max;
        mean_out{iC,iF} = binned_mean;
    end
end
fig = figure('WindowState','maximized');
for iF = 1:length(f_list)
    
    % Plot mean
    subplot(2,length(f_list),iF)
    hold on;
    for iC = 1:length(max_out)
        plot(mean_out{iC,iF}, 'color',c_list{iF}, 'LineWidth', 0.25);
    end
    temp_ppc = cell2mat(mean_out(:,iF));
%     temp_ppc(~isfinite(temp_ppc)) = NaN;     % Convert Inf and -Inf to nan
    plot(nanmean(temp_ppc), 'color', 'black', 'LineWidth', 3);
    xticks([1.5, 2.5, 3.5, 4.5, 5.5])
    xticklabels({'-2', '-1', '0', '1', '2'})
    xlim([0.85 7.15])
    ylabel(' PPC')
    xlabel('Z-scored FR')
    title(sprintf("Mean: %d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})


    % Plot max
    subplot(2,length(f_list),iF+length(f_list))
    hold on;
    for iC = 1:length(max_out)
        plot(max_out{iC,iF}, 'color',c_list{iF}, 'LineWidth', 0.25);
    end
    temp_ppc = cell2mat(max_out(:,iF));
%     temp_ppc(~isfinite(temp_ppc)) = NaN;     % Convert Inf and -Inf to nan
    plot(nanmean(temp_ppc), 'color', 'black', 'LineWidth', 3);
    xticks([1.5, 2.5, 3.5, 4.5, 5.5])
    xticklabels({'-2', '-1', '0', '1', '2'})
    xlim([0.85 7.15])
    ylabel('PPC')
    xlabel('Z-scored FR')
    title(sprintf("Max: %d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
end



