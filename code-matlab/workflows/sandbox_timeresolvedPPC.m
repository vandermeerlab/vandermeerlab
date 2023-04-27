%% Script to calculate timeresolved PPC for cells that have significant phase locking

% setup parameters
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items

purple = [0.7 0 1]; % Color code for purple
fsi_summary = table("label","% trials with Feeder2", ...
    "delta_v_corr","delta_d_corr", "c_delta_v_corr","c_delta_d_corr", ...
    "theta_v_corr","theta_d_corr", "c_theta_v_corr","c_theta_d_corr", ...
    "beta_v_corr","beta_d_corr", "c_beta_v_corr","c_beta_d_corr", ...
    "lowgamma_v_corr","lowgamma_d_corr", "c_lowgamma_v_corr","c_lowgamma_d_corr", ...
    "highgamma_v_corr","highgamma_d_corr", "c_highgamma_v_corr","c_highgamma_d_corr");
msn_summary = table("label","% trials with Feeder2", ...
    "delta_v_corr","delta_d_corr", "c_delta_v_corr","c_delta_d_corr", ...
    "theta_v_corr","theta_d_corr", "c_theta_v_corr","c_theta_d_corr", ...
    "beta_v_corr","beta_d_corr", "c_beta_v_corr","c_beta_d_corr", ...
    "lowgamma_v_corr","lowgamma_d_corr", "c_lowgamma_v_corr","c_lowgamma_d_corr", ...
    "highgamma_v_corr","highgamma_d_corr", "c_highgamma_v_corr","c_highgamma_d_corr");

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
    
    % Taking a longer window than what we need (+/-5 seconds) to avoid edge
    % effects
    w_start = rt1 - 7;
    w_end =  rt1 + 7;
    % If the last trial time exceeds end of experiment, discard that trial
    if w_end(end) < ExpKeys.TimeOffTrack
        w_end(end) = []; w_start(end) = [];
    end
    rt_iv = iv(w_start, w_end);

    % Save a list of the rt2 that are within the trial-window
    rt2 = getRewardTimes2();
    valid_rt2 = [];
    for iT = 1:length(rt2)
        vidx = find((rt2(iT) >= (w_start+2)) & (rt2(iT) <= (w_end-2)));
%         disp(vidx)
        assert(length(vidx)<2,sprintf("You have got problems with trial defintions in %s", label));
        if ~isempty(vidx)
            valid_rt2 = [valid_rt2, rt2(iT) - w_start(vidx) - 2]; % the -2 is to account for the hack we are doing
        end  
    end

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
        warning('Cell %s has nans in its STS',sd.S.label{1});
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

    % Calculate time-resolved PPC
    cfg_ppc2 = cfg_ppc;
    cfg_ppc2.latency = 'maxperiod';
    cfg_ppc2.timwin = 2;
    cfg_ppc2.winstepsize = 0.1;
    this_tr_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc2, this_sts);
    mean_tr_fr = squeeze(this_tr_ppc.nspikes(:,1,:))'/size(near_data.trial,2);
    % running dxdt before subsetting to avoid edge effects
    dfr = dxdt(this_tr_ppc.time, mean_tr_fr, 'postSmoothing', 0);
    % dxdt doesn't work on 2D arrays, so need to run it for each frequency band
    this_tr_ppc.f_ppc = zeros(length(f_list), length(this_tr_ppc.time));
    this_tr_ppc.d_ppc = zeros(length(f_list), length(this_tr_ppc.time));
    temp_ppc = squeeze(this_tr_ppc.ppc0);
    for iB = 1:length(f_list)
        f_idx = find(round(this_freq) >= f_list{iB}(1) & ...
            round(this_freq) <= f_list{iB}(2));
        this_tr_ppc.f_ppc(iB,:) = mean(temp_ppc(f_idx,:));
        this_tr_ppc.d_ppc(iB,:) = dxdt(this_tr_ppc.time, this_tr_ppc.f_ppc(iB,:), 'postSmoothing', 0);
    end

    % subset the results to only keep the trial +/- 5 sec window
    keep = ((this_tr_ppc.time - 7) >= -5) & ((this_tr_ppc.time - 7) <= 5);
    this_tr_ppc.time = this_tr_ppc.time(keep) - 2;
    this_tr_ppc.ppc0 = this_tr_ppc.ppc0(:,:,keep);
    this_tr_ppc.nspikes = this_tr_ppc.nspikes(:,:,keep);
    this_tr_ppc.f_ppc = this_tr_ppc.f_ppc(:,keep);
    this_tr_ppc.d_ppc = this_tr_ppc.d_ppc(:,keep);
    mean_tr_fr = mean_tr_fr(keep);
    dfr = dfr(keep);
    
    % Find the earliest reward2 firing time point and save clean results
    first_rt2 = nearest_idx3(valid_rt2, this_tr_ppc.time);
    first_rt2 = min(first_rt2);
    clean_tr_fr = mean_tr_fr;
    clean_tr_ppc = this_tr_ppc;
    clean_dfr = dfr;
    if ~isempty(first_rt2)
        clean_tr_fr = clean_tr_fr(1:first_rt2);
        clean_tr_ppc.time = clean_tr_ppc.time(1:first_rt2);
        clean_tr_ppc.ppc0 = clean_tr_ppc.ppc0(:,:,1:first_rt2);
        clean_tr_ppc.nspikes = clean_tr_ppc.nspikes(:,:,1:first_rt2);
        clean_dfr = clean_dfr(1:first_rt2);
        clean_tr_ppc.f_ppc = clean_tr_ppc.f_ppc(:,1:first_rt2);
        clean_tr_ppc.d_ppc = clean_tr_ppc.d_ppc(:,1:first_rt2);
    end

    % Plotting
    this_bands = find(this_sig);
    tp = length(this_bands) + 1;
    fig = figure('WindowState', 'maximized');
    % Plot whole window 
    ax = subplot(tp, 2, 1);
    hold on
    if ~isempty(valid_rt2)
        xline(valid_rt2 - 5, '--red')
    end
    plot(this_tr_ppc.time-5, mean_tr_fr, 'Color', purple)
    yyaxis(ax, 'right');
    plot(this_tr_ppc.time-5, dfr, 'black')
    ax.XAxis.Label.String = 'Time';
    ax.XAxis.Label.FontSize = 12;
    ax.YAxis(1).Label.String = 'Inst. F.R';
    ax.YAxis(1).Label.FontSize = 12;
    ax.YAxis(1).Color = purple;
    ax.YAxis(2).Label.Interpreter = 'tex';
    ax.YAxis(2).Label.String = '\delta F.R';
    ax.YAxis(2).Label.FontSize = 12;
    ax.YAxis(2).Color = 'black';
    pct_rt2 = (100*length(valid_rt2))/size(near_data.trial,2);
    this_row = {label, pct_rt2, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, ...
        nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
    ax.Title.String = sprintf('%.2f%% trials with Feeder2 firing',pct_rt2);
    ax.Title.Color = 'red';
    ax.Title.Interpreter = 'none';
    
    % Plot clean window
    ax = subplot(tp,2, 2);
    hold on;
    plot(clean_tr_ppc.time-5, clean_tr_fr, 'Color', purple)
    yyaxis(ax, 'right')
    plot(clean_tr_ppc.time-5, clean_dfr, 'black')
    ax.XAxis.Label.String = 'Time';
    ax.XAxis.Label.FontSize = 12;
    ax.XLim = [-5,round(clean_tr_ppc.time(end),1)-5];
    ax.YAxis(1).Label.String = 'Inst. F.R';
    ax.YAxis(1).Label.FontSize = 12;
    ax.YAxis(1).Color = purple;
    ax.YAxis(2).Label.Interpreter = 'tex';
    ax.YAxis(2).Label.String = '\delta F.R';
    ax.YAxis(2).Label.FontSize = 12;
    ax.YAxis(2).Color = 'black';
    ax.Title.Interpreter = 'none';
    ax.Title.String = 'No Reward2 firing';

    for iB = 1:length(this_bands)
        % Plot whole window
        f_ppc = this_tr_ppc.f_ppc(this_bands(iB),:);
        ax = subplot(tp,2,(2*iB)+1);
        hold on
        plot(this_tr_ppc.time-5, f_ppc, 'Color', c_list{this_bands(iB)})
        ax.XAxis.Label.String = 'Time';
        ax.XAxis.Label.FontSize = 12;
        v_corr = corrcoef(f_ppc, mean_tr_fr);
        v_corr = v_corr(1,2); % Correlation between PPC and FR

        yyaxis(ax, 'right');
        d_ppc = this_tr_ppc.d_ppc(this_bands(iB),:);
        plot(this_tr_ppc.time-5, d_ppc, 'black');
        ax.YAxis(1).Label.String = 'PPC';
        ax.YAxis(1).Label.FontSize = 12;
        ax.YAxis(1).Exponent = 0;
        ax.YAxis(1).Color = c_list{this_bands(iB)};
        ax.YAxis(2).Label.Interpreter = 'tex';
        ax.YAxis(2).Label.String = '\delta PPC';
        ax.YAxis(2).Label.FontSize = 12;
        ax.YAxis(2).Exponent = 0;
        ax.YAxis(2).Color = 'black';
        d_corr = corrcoef(d_ppc, dfr);
        d_corr = d_corr(1,2); % Correlation between delta PPC and delta FR
        ax.Title.Interpreter = 'tex';
        ax.Title.Color = c_list{this_bands(iB)};
        ax.Title.String = sprintf("%d Hz - %d Hz, v-corr:%.2f, {\\color{black}d-corr:%.2f}", ...
            f_list{this_bands(iB)}(1), f_list{this_bands(iB)}(2), v_corr, d_corr);
        
        % Plot clean window
        clean_f_ppc = clean_tr_ppc.f_ppc(this_bands(iB),:);
        ax = subplot(tp,2,(2*iB)+2);
        hold on
        plot(clean_tr_ppc.time-5, clean_f_ppc, 'Color', c_list{this_bands(iB)})
        ax.XAxis.Label.String = 'Time';
        ax.XAxis.Label.FontSize = 12;
        ax.XLim = [-5,round(clean_tr_ppc.time(end),1)-5];
        clean_v_corr = corrcoef(clean_f_ppc, clean_tr_fr);
        clean_v_corr = clean_v_corr(1,2); % Correlation between PPC and FR

        yyaxis(ax, 'right');
        clean_d_ppc = clean_tr_ppc.d_ppc(this_bands(iB),:);
        plot(clean_tr_ppc.time-5, clean_d_ppc, 'black');
        ax.YAxis(1).Label.String = 'PPC';
        ax.YAxis(1).Label.FontSize = 12;
        ax.YAxis(1).Exponent = 0;
        ax.YAxis(1).Color = c_list{this_bands(iB)};
        ax.YAxis(2).Label.Interpreter = 'tex';
        ax.YAxis(2).Label.String = '\delta PPC';
        ax.YAxis(2).Label.FontSize = 12;
        ax.YAxis(2).Exponent = 0;
        ax.YAxis(2).Color = 'black';
        clean_d_corr = corrcoef(clean_d_ppc, clean_dfr);
        clean_d_corr = clean_d_corr(1,2); % Correlation between PPC and FR
        ax.Title.Interpreter = 'tex';
        ax.Title.Color = c_list{this_bands(iB)};
        ax.Title.String = sprintf("%d Hz - %d Hz, v-corr:%.2f, {\\color{black}d-corr:%.2f}", ...
            f_list{this_bands(iB)}(1), f_list{this_bands(iB)}(2), clean_v_corr, clean_d_corr);
        
        % Saving results
        this_row{this_bands(iB)*4-1} = v_corr;
        this_row{this_bands(iB)*4} = d_corr;
        this_row{this_bands(iB)*4+1} = clean_v_corr;
        this_row{this_bands(iB)*4+2} = clean_d_corr;
    end
    fsi_summary = [fsi_summary; this_row];
    sgtitle(sprintf('FSI: %s', label), 'Interpreter', 'none')
    WriteFig(fig, strcat('D:\RandomVstrAnalysis\PPCbyFR\FSI_',label,'_FRbyPPC'),1)
    close;
end
writetable(fsi_summary, strcat('D:\RandomVstrAnalysis\PPCbyFR\PhaseLocked\','fsi_summary.xls'));
%% Cycle through MSNs next
for iC = 1:size(sig3_msn,1)
    this_sig = sig3_msn{iC,2};
    % Skip this cell if it is not significantly phase locked to any
    % frequency band
    if ~any(this_sig)
        continue
    end
    label = sig3_msn{iC,1};
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
    
    % Taking a longer window than what we need (+/-5 seconds) to avoid edge
    % effects
    w_start = rt1 - 7;
    w_end =  rt1 + 7;
    % If the last trial time exceeds end of experiment, discard that trial
    if w_end(end) < ExpKeys.TimeOffTrack
        w_end(end) = []; w_start(end) = [];
    end
    rt_iv = iv(w_start, w_end);

    % Save a list of the rt2 that are within the trial-window
    rt2 = getRewardTimes2();
    valid_rt2 = [];
    for iT = 1:length(rt2)
        vidx = find((rt2(iT) >= (w_start+2)) & (rt2(iT) <= (w_end-2)));
%         disp(vidx)
        assert(length(vidx)<2,sprintf("You have got problems with trial defintions in %s", label));
        if ~isempty(vidx)
            valid_rt2 = [valid_rt2, rt2(iT) - w_start(vidx) - 2]; % the -2 is to account for the hack we are doing
        end  
    end

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
        warning('Cell %s has nans in its STS',sd.S.label{1});
    end
    this_freq = this_sts.freq;
    this_sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
    flag_nansts = this_flag;
    
    % Calculate and save trialwise_ppc and use subsampled measures
    % for near Reward Trials. However, if an MSN has fewer spikes
    % than at least one co-recorded MSN, do NOT
    % subsample and include all the spikes
    
    % Calculate non time-resolved PPC
    cfg_ppc               = [];
    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
    cfg_ppc.spikechannel  = this_sts.label;
    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
    cfg_ppc.avgoverchan   = 'weighted';
    this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc, this_sts);

    % Calculate time-resolved PPC
    cfg_ppc2 = cfg_ppc;
    cfg_ppc2.latency = 'maxperiod';
    cfg_ppc2.timwin = 2;
    cfg_ppc2.winstepsize = 0.1;
    this_tr_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc2, this_sts);
    mean_tr_fr = squeeze(this_tr_ppc.nspikes(:,1,:))'/size(near_data.trial,2);
    % running dxdt before subsetting to avoid edge effects
    dfr = dxdt(this_tr_ppc.time, mean_tr_fr, 'postSmoothing', 0);
    % dxdt doesn't work on 2D arrays, so need to run it for each frequency band
    this_tr_ppc.f_ppc = zeros(length(f_list), length(this_tr_ppc.time));
    this_tr_ppc.d_ppc = zeros(length(f_list), length(this_tr_ppc.time));
    temp_ppc = squeeze(this_tr_ppc.ppc0);
    for iB = 1:length(f_list)
        f_idx = find(round(this_freq) >= f_list{iB}(1) & ...
            round(this_freq) <= f_list{iB}(2));
        this_tr_ppc.f_ppc(iB,:) = mean(temp_ppc(f_idx,:));
        this_tr_ppc.d_ppc(iB,:) = dxdt(this_tr_ppc.time, this_tr_ppc.f_ppc(iB,:), 'postSmoothing', 0);
    end

    % subset the results to only keep the trial +/- 5 sec window
    keep = ((this_tr_ppc.time - 7) >= -5) & ((this_tr_ppc.time - 7) <= 5);
    this_tr_ppc.time = this_tr_ppc.time(keep) - 2;
    this_tr_ppc.ppc0 = this_tr_ppc.ppc0(:,:,keep);
    this_tr_ppc.nspikes = this_tr_ppc.nspikes(:,:,keep);
    this_tr_ppc.f_ppc = this_tr_ppc.f_ppc(:,keep);
    this_tr_ppc.d_ppc = this_tr_ppc.d_ppc(:,keep);
    mean_tr_fr = mean_tr_fr(keep);
    dfr = dfr(keep);
    
    % Find the earliest reward2 firing time point and save clean results
    first_rt2 = nearest_idx3(valid_rt2, this_tr_ppc.time);
    first_rt2 = min(first_rt2);
    clean_tr_fr = mean_tr_fr;
    clean_tr_ppc = this_tr_ppc;
    clean_dfr = dfr;
    if ~isempty(first_rt2)
        clean_tr_fr = clean_tr_fr(1:first_rt2);
        clean_tr_ppc.time = clean_tr_ppc.time(1:first_rt2);
        clean_tr_ppc.ppc0 = clean_tr_ppc.ppc0(:,:,1:first_rt2);
        clean_tr_ppc.nspikes = clean_tr_ppc.nspikes(:,:,1:first_rt2);
        clean_dfr = clean_dfr(1:first_rt2);
        clean_tr_ppc.f_ppc = clean_tr_ppc.f_ppc(:,1:first_rt2);
        clean_tr_ppc.d_ppc = clean_tr_ppc.d_ppc(:,1:first_rt2);
    end

    % Plotting
    this_bands = find(this_sig);
    tp = length(this_bands) + 1;
    fig = figure('WindowState', 'maximized');
    % Plot whole window 
    ax = subplot(tp, 2, 1);
    hold on
    if ~isempty(valid_rt2)
        xline(valid_rt2 - 5, '--red')
    end
    plot(this_tr_ppc.time-5, mean_tr_fr, 'Color', purple)
    yyaxis(ax, 'right');
    plot(this_tr_ppc.time-5, dfr, 'black')
    ax.XAxis.Label.String = 'Time';
    ax.XAxis.Label.FontSize = 12;
    ax.YAxis(1).Label.String = 'Inst. F.R';
    ax.YAxis(1).Label.FontSize = 12;
    ax.YAxis(1).Color = purple;
    ax.YAxis(2).Label.Interpreter = 'tex';
    ax.YAxis(2).Label.String = '\delta F.R';
    ax.YAxis(2).Label.FontSize = 12;
    ax.YAxis(2).Color = 'black';
    pct_rt2 = (100*length(valid_rt2))/size(near_data.trial,2);
    this_row = {label, pct_rt2, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, ...
        nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
    ax.Title.String = sprintf('%.2f%% trials with Feeder2 firing',pct_rt2);
    ax.Title.Color = 'red';
    ax.Title.Interpreter = 'none';
    
    % Plot clean window
    ax = subplot(tp,2, 2);
    hold on;
    plot(clean_tr_ppc.time-5, clean_tr_fr, 'Color', purple)
    yyaxis(ax, 'right')
    plot(clean_tr_ppc.time-5, clean_dfr, 'black')
    ax.XAxis.Label.String = 'Time';
    ax.XAxis.Label.FontSize = 12;
    ax.XLim = [-5,round(clean_tr_ppc.time(end),1)-5];
    ax.YAxis(1).Label.String = 'Inst. F.R';
    ax.YAxis(1).Label.FontSize = 12;
    ax.YAxis(1).Color = purple;
    ax.YAxis(2).Label.Interpreter = 'tex';
    ax.YAxis(2).Label.String = '\delta F.R';
    ax.YAxis(2).Label.FontSize = 12;
    ax.YAxis(2).Color = 'black';
    ax.Title.Interpreter = 'none';
    ax.Title.String = 'No Reward2 firing';

    for iB = 1:length(this_bands)
        % Plot whole window
        f_ppc = this_tr_ppc.f_ppc(this_bands(iB),:);
        ax = subplot(tp,2,(2*iB)+1);
        hold on
        plot(this_tr_ppc.time-5, f_ppc, 'Color', c_list{this_bands(iB)})
        ax.XAxis.Label.String = 'Time';
        ax.XAxis.Label.FontSize = 12;
        v_corr = corrcoef(f_ppc, mean_tr_fr);
        v_corr = v_corr(1,2); % Correlation between PPC and FR

        yyaxis(ax, 'right');
        d_ppc = this_tr_ppc.d_ppc(this_bands(iB),:);
        plot(this_tr_ppc.time-5, d_ppc, 'black');
        ax.YAxis(1).Label.String = 'PPC';
        ax.YAxis(1).Label.FontSize = 12;
        ax.YAxis(1).Exponent = 0;
        ax.YAxis(1).Color = c_list{this_bands(iB)};
        ax.YAxis(2).Label.Interpreter = 'tex';
        ax.YAxis(2).Label.String = '\delta PPC';
        ax.YAxis(2).Label.FontSize = 12;
        ax.YAxis(2).Exponent = 0;
        ax.YAxis(2).Color = 'black';
        d_corr = corrcoef(d_ppc, dfr);
        d_corr = d_corr(1,2); % Correlation between delta PPC and delta FR
        ax.Title.Interpreter = 'tex';
        ax.Title.Color = c_list{this_bands(iB)};
        ax.Title.String = sprintf("%d Hz - %d Hz, v-corr:%.2f, {\\color{black}d-corr:%.2f}", ...
            f_list{this_bands(iB)}(1), f_list{this_bands(iB)}(2), v_corr, d_corr);
        
        % Plot clean window
        clean_f_ppc = clean_tr_ppc.f_ppc(this_bands(iB),:);
        ax = subplot(tp,2,(2*iB)+2);
        hold on
        plot(clean_tr_ppc.time-5, clean_f_ppc, 'Color', c_list{this_bands(iB)})
        ax.XAxis.Label.String = 'Time';
        ax.XAxis.Label.FontSize = 12;
        ax.XLim = [-5,round(clean_tr_ppc.time(end),1)-5];
        clean_v_corr = corrcoef(clean_f_ppc, clean_tr_fr);
        clean_v_corr = clean_v_corr(1,2); % Correlation between PPC and FR

        yyaxis(ax, 'right');
        clean_d_ppc = clean_tr_ppc.d_ppc(this_bands(iB),:);
        plot(clean_tr_ppc.time-5, clean_d_ppc, 'black');
        ax.YAxis(1).Label.String = 'PPC';
        ax.YAxis(1).Label.FontSize = 12;
        ax.YAxis(1).Exponent = 0;
        ax.YAxis(1).Color = c_list{this_bands(iB)};
        ax.YAxis(2).Label.Interpreter = 'tex';
        ax.YAxis(2).Label.String = '\delta PPC';
        ax.YAxis(2).Label.FontSize = 12;
        ax.YAxis(2).Exponent = 0;
        ax.YAxis(2).Color = 'black';
        clean_d_corr = corrcoef(clean_d_ppc, clean_dfr);
        clean_d_corr = clean_d_corr(1,2); % Correlation between PPC and FR
        ax.Title.Interpreter = 'tex';
        ax.Title.Color = c_list{this_bands(iB)};
        ax.Title.String = sprintf("%d Hz - %d Hz, v-corr:%.2f, {\\color{black}d-corr:%.2f}", ...
            f_list{this_bands(iB)}(1), f_list{this_bands(iB)}(2), clean_v_corr, clean_d_corr);
        
        % Saving results
        this_row{this_bands(iB)*4-1} = v_corr;
        this_row{this_bands(iB)*4} = d_corr;
        this_row{this_bands(iB)*4+1} = clean_v_corr;
        this_row{this_bands(iB)*4+2} = clean_d_corr;
    end
    msn_summary = [msn_summary; this_row];
    sgtitle(sprintf('MSN: %s', label), 'Interpreter', 'none')
    WriteFig(fig, strcat('D:\RandomVstrAnalysis\PPCbyFR\MSN_',label,'_FRbyPPC'),1)
    close;
end
writetable(msn_summary, strcat('D:\RandomVstrAnalysis\PPCbyFR\PhaseLocked\','msn_summary.xls'));
%%
 

