%% script to generate trial wise PPC for MSNs ans FSIs
% The spectra for the FSIs are calculated after multiple rounds of
% subsampling
%% setup
clear;
% Setup rng seed for reproducibility
rng(4994);
cd('E:\ADRLabData');
% cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData');
please = [];
please.rats = {'R117','R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.output_dir = 'D:\RandomVstrAnalysis\temp_no_thresh';
% cfg_in.output_dir = '/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp';
cfg_in.incl_types = [1, 2];
cfg_in.nMinSpikes1 = 400; % For on track
cfg_in.nMinSpikes2 = 400; % For near
cfg_in.nMinSpikes3 = 150; % For lfr, hfr, p1 and p2
cfg_in.nControlSplits = 100;
cfg_in.num_subsamples = 1000;
cfg_in.minTrialSpikes = 2; % Not thresholding


%%
% Top level loop which calls the main function for all the sessions
for iS = 1:length(cfg_in.fd) % for each session...
    cfg_in.iS = iS;
    pushdir(cfg_in.fd{iS});
    generateSTS(cfg_in); % do the business
    popdir;
end % of sessions

%%
% Main function to generate spike_triggered_spectra
function od = generateSTS(cfg_in)

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
        
    lfp_tt = regexp(cfg.fc, 'CSC\d+', 'match');
    lfp_tt = str2double(lfp_tt{1}{1}(4:end)); % need this to skip cells from same tt (could make into function)
    fprintf('LFP ttno is %d\n', lfp_tt);
    
    % params
    cfg_master = []; % overall params
    cfg_master.dt = 0.001;
    cfg_master.ccMethod = 'MvdM'; % cell type classification method
    cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
    cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
    cfg_master.iS = 1; % current session number out of fd list, get this from input cfg
    cfg_master.fd = []; % full list of session fd's, get this from input cfg
    cfg_master.fd_extra = []; % get this from input cfg
    cfg_master.write_output = 0;
    cfg_master.output_prefix = 'sts_';
    cfg_master.output_dir = 'C:\temp';
    cfg_master.incl_cell_types = [1,2]; %cell types to be included
    cfg_master.num_subsamples = 1000;
    cfg_master = ProcessConfig(cfg_master,cfg_in);
    
    % spikes
    sd.S = LoadSpikesTarget(cfg_master);
    
    % Restrict spikes to only OnTrack
    sd.S = restrict(sd.S, iv(ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack));
    
    % Categorize cells and add tetrode depths
    cfg_wv = []; cfg_wv.cMethod = cfg_master.ccMethod;
    s_out = CategorizeStriatumWave(cfg_wv, sd.S);

    s_out.unit = [s_out.other s_out.msn s_out.fsi];
    s_out.ident = [zeros(1, length(s_out.other)) ones(1, length(s_out.msn)) repmat(2, 1, length(s_out.fsi))];

    cfg_tt = []; cfg_tt.verbose = 1;
    cfg_tt.this_rat_ID = cfg_master.fd_extra.ratID_num(cfg_master.iS);
    cfg_tt.this_date = cfg_master.fd_extra.fd_date_num(cfg_master.iS);

    for iC = 1:length(sd.S.t)
        % Read the spike files into field trip format
        sd.S.ft_spikes(iC) = ft_read_spike(sd.S.label{iC});
        sd.S.usr.cell_type(iC) = s_out.ident((s_out.unit == iC));
        sd.S.usr.tetrodeDepths(iC) = ExpKeys.TetrodeDepths(sd.S.usr.tt_num(iC));
        cfg_tt.ttno = sd.S.usr.tt_num(iC);
        [sd.S.usr.distanceTurned(iC), prev_fd] = DistanceTurned(cfg_tt, cfg_master.fd, cfg_master.fd_extra);
        cfg_tt.verbose = 0;
    end
     
    % correlate with previous session waveforms if available
    if isempty(prev_fd) % no previous day available
        sd.S.usr.duplicate = zeros(size(sd.S.usr.tt_num));
    else
        pushdir(prev_fd);
        S2 = LoadSpikes([]);
        nSpikes = cellfun(@length, S2.t); keep = nSpikes >= cfg_master.nMinSpikes1;
        S2 = SelectTS([], S2, keep);

        s_out2 = CategorizeStriatumWave(cfg_wv, S2);
        s_out = CalcWVDistances([], s_out, s_out2); % add comparison with previous day's waveforms

        popdir; 

        % for each cell in current session, decide if duplicate
        for iC = 1:length(sd.S.t)

            this_tt_no = sd.S.usr.tt_num(iC);
            prev_day_cells = find(S2.usr.tt_num == this_tt_no);

            if isempty(prev_day_cells) % no cells recorded fron this tt in previous session
                sd.S.usr.duplicate(iC) = 0;
            else % previous day cells found
                temp_corr = s_out.corr(iC, prev_day_cells);
                temp_peakn = s_out.peakdiffn(iC, prev_day_cells);

                if temp_corr > cfg_master.maxPrevCorr & abs(temp_peakn) < cfg_master.maxPeakn % wv correlation big, peak difference small
                    sd.S.usr.duplicate(iC) = 1;
                else
                    sd.S.usr.duplicate(iC) = 0;
                end
            end
        end
    end % of previous day available checks
    
    % PLEASE SEE: Use the sd.S.usr.duplicate field if you are calculating
    % it anyway!!
    % Keep only non duplicate cells and those which were read by FieldTrip 
    % correctly
    % Also Get rid of cells of exc_types or spikes < nMinSpikes1;
    % Also get rid of cells on the LFP tt
    keep = false(1,length(sd.S.t));
    for iK = 1:length(keep)
        for iT = 1:length(cfg_master.incl_cell_types)
           keep(iK) = keep(iK) | (cfg_master.incl_cell_types(iT) == sd.S.usr.cell_type(iK)); 
        end
        keep(iK) = keep(iK) & (length(sd.S.t{iK}) > cfg_master.nMinSpikes1) & ~(sd.S.usr.tt_num(iK) == lfp_tt);
    end
    keep = keep & ~sd.S.usr.duplicate;% & sd.S.ft_spk_valid;
    sd.S = SelectTS([], sd.S, keep);
    sd.S.ft_spikes = sd.S.ft_spikes(keep);
    od.cell_type = sd.S.usr.cell_type;
    od.tt_id = sd.S.usr.tt_num;
    od.label = sd.S.label;
   
    % Calculate spectral measures for all MSNs
    all_msn = find(od.cell_type == 1);
    for iM = 1:length(all_msn)
        iC  = all_msn(iM);
        % Calculate and save STA
        cfg_ft.timwin = [-0.5 0.5];
        cfg_ft.spikechannel = sd.S.ft_spikes(iC).label{1};
        cfg_ft.channel = ft_csc.label(1);
        this_data = ft_appendspike([], ft_csc, sd.S.ft_spikes(iC));
 
        % Block of code to divide recordings session into near trials

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
        spk_count1 = length(restrict(sd.S, rt_iv).t{iC});
        spk_count2 = 0;
        % Save trial_wise spike conunt
        near_trialwise_spk_count = zeros(1, length(near_data.trial));
        for iT = 1:length(near_data.trial)
           near_trialwise_spk_count(iT) = sum(near_data.trial{iT}(2,:));
           spk_count2 = spk_count2 + near_trialwise_spk_count(iT);
        end
        % Skip if no spikes present!
        if spk_count2 <  cfg_master.nMinSpikes2
            od.msn_res.near_spec{iM}.flag_tooFewSpikes = true;
        else
            od.msn_res.near_spec{iM}.spk_count = spk_count2;
            od.msn_res.near_spec{iM}.flag_tooFewSpikes = false;
            if spk_count1 ~= spk_count2
                this_flag = true;
                warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in near-Reward trials', ...
                    spk_count1, spk_count2)

            end
            od.msn_res.near_spec{iM}.trialwise_spk_count = near_trialwise_spk_count;
            od.msn_res.near_spec{iM}.flag_unequalSpikes = this_flag;

            % Extract All Spike IDs
            trial_wise_spike = cell(1, length(near_data.trial));
            last_spk_ct = 0;
            for iT = 1:length(near_data.trial)
                this_spk_ct = near_trialwise_spk_count(iT);
                trial_wise_spike{iT} = last_spk_ct + 1 : last_spk_ct + this_spk_ct;
                last_spk_ct = last_spk_ct + this_spk_ct;
            end

            % Calculate and save STS for all the trials
            cfg_sts.method = 'mtmconvol';
            cfg_sts.foi = 1:1:100;
            cfg_sts.t_ftimwin = 5./cfg_sts.foi;
            cfg_sts.taper = 'hanning';
            cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
            cfg_sts.channel = near_data.label{1};
            this_sts = ft_spiketriggeredspectrum(cfg_sts, near_data);
            this_flag = false;
            % Display warning to show that there were Nans in this calculation
            if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                this_flag = true;
                warning('Cell %s has nans in its STS',sd.S.label{iC});
            end
            od.msn_res.near_spec{iM}.freqs = this_sts.freq;
            od.msn_res.near_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
            od.msn_res.near_spec{iM}.flag_nansts = this_flag;

            % Calculate and save trialwise_ppc
            near_trialwise_ppc = zeros(length(near_data.trial), length(od.msn_res.near_spec{iM}.freqs));
            cfg_ppc               = [];
            cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
            cfg_ppc.spikechannel  = this_sts.label;
            cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
            cfg_ppc.avgoverchan   = 'weighted';
            this_flag = false;
            % If a trial has only one spike, ppc can't be calculated!
            % Need at at least 2 spikes!
            for iT = 1:length(near_data.trial)
                 if od.msn_res.near_spec{iM}.trialwise_spk_count(iT) < cfg_in.minTrialSpikes
                     continue;
                 else
                    this_trial_sts = this_sts;
                    this_trial_spks = trial_wise_spike{iT};
                    this_trial_sts.fourierspctrm{1} = this_trial_sts.fourierspctrm{1}(this_trial_spks,:,:);
                    this_trial_sts.time{1} = this_trial_sts.time{1}(this_trial_spks,:);
                    this_trial_sts.trial{1} = this_trial_sts.trial{1}(this_trial_spks,:);
                    this_trial_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc,this_trial_sts);
                    if ~isempty(find(isnan(this_trial_ppc.ppc0),1))
                        this_flag = true;
                        break;
                    end
                    near_trialwise_ppc(iT,:) = this_trial_ppc.ppc0;
                 end
            end
            od.msn_res.near_spec{iM}.flag_nanppc = this_flag;
            if ~this_flag
                od.msn_res.near_spec{iM}.trial_wise_ppc = near_trialwise_ppc;
            end
        end 
    end
    
    % Get the msn distributions
    od.msn_near_dist = [];
    
    % Subsample only if all the splits are problem free
    % Also convert all numbers less than 2 to 0 because you need at least 2
    % spikes in a trial to calculate PPC
    if isfield(od,'msn_res') && isfield(od.msn_res, 'near_spec')
        for iM = 1:length(od.msn_res.near_spec)
            if ~od.msn_res.near_spec{iM}.flag_tooFewSpikes & ...
                ~od.msn_res.near_spec{iM}.flag_nanppc
                    temp_spk_count = od.msn_res.near_spec{iM}.trialwise_spk_count';
                    temp_spk_count(temp_spk_count < 2) = 0;
                    od.msn_near_dist = [od.msn_near_dist, temp_spk_count];
            end
        end
    end
    
    % Calculate spectral measures for all FSIs
    all_fsi = find(od.cell_type == 2);
    for iM = 1:length(all_fsi)   
        iC  = all_fsi(iM);
        if isempty(od.msn_near_dist)
            continue;
        end
        % Calculate and save STA
        cfg_ft.timwin = [-0.5 0.5];
        cfg_ft.spikechannel = sd.S.ft_spikes(iC).label{1};
        cfg_ft.channel = ft_csc.label(1);
        this_data = ft_appendspike([], ft_csc, sd.S.ft_spikes(iC));
       
        % Block of code to divide recordings session into near and away trials

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
        spk_count1 = length(restrict(sd.S, rt_iv).t{iC});
        spk_count2 = 0;
        % Save trial_wise spike conunt
        near_trialwise_spk_count = zeros(1, length(near_data.trial));
        for iT = 1:length(near_data.trial)
           near_trialwise_spk_count(iT) = sum(near_data.trial{iT}(2,:));
           spk_count2 = spk_count2 + near_trialwise_spk_count(iT);
        end
        % Skip if no spikes present!
        if spk_count2 <  cfg_master.nMinSpikes2
            od.fsi_res.near_spec{iM}.flag_tooFewSpikes = true;
        else
            od.fsi_res.near_spec{iM}.spk_count = spk_count2;
            od.fsi_res.near_spec{iM}.flag_tooFewSpikes = false;
            if spk_count1 ~= spk_count2
                this_flag = true;
                warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in near-Reward trials', ...
                    spk_count1, spk_count2)
            end
            od.fsi_res.near_spec{iM}.trialwise_spk_count = near_trialwise_spk_count;
            od.fsi_res.near_spec{iM}.flag_unequalSpikes = this_flag;
            
            % Extract All Spike IDs
            trial_wise_spike = cell(1, length(near_data.trial));
            last_spk_ct = 0;
            for iT = 1:length(near_data.trial)
                this_spk_ct = near_trialwise_spk_count(iT);
                trial_wise_spike{iT} = last_spk_ct + 1 : last_spk_ct + this_spk_ct;
                last_spk_ct = last_spk_ct + this_spk_ct;
            end

            % Calculate and save STS for all spikes
            cfg_sts.method = 'mtmconvol';
            cfg_sts.foi = 1:1:100;
            cfg_sts.t_ftimwin = 5./cfg_sts.foi;
            cfg_sts.taper = 'hanning';
            cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
            cfg_sts.channel = near_data.label{1};
            this_sts = ft_spiketriggeredspectrum(cfg_sts, near_data);
            this_flag = false;
            % Display warning to show that there were Nans in this calculation
            if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                this_flag = true;
                warning('Cell %s has nans in its STS',sd.S.label{iC});
            end
            od.fsi_res.near_spec{iM}.freqs = this_sts.freq;
            od.fsi_res.near_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
            od.fsi_res.near_spec{iM}.flag_nansts = this_flag;
            
            % Calculate and save trialwise_ppc and use subsampled measures
            % for near Reward Trials. However, if an FSI has fewer spikes
            % than at least one co-recorded MSN, do NOT
            % subsample and include all the spikes
            
            % Generate num_subsamples versions of trial_wise_ppc
            od.fsi_res.near_spec{iM}.trialwise_ppc = nan(cfg_master.num_subsamples, ...
                length(near_data.trial), length(od.fsi_res.near_spec{iM}.freqs));
            cfg_ppc               = [];
            cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
            cfg_ppc.spikechannel  = this_sts.label;
            cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
            cfg_ppc.avgoverchan   = 'weighted';
            this_flag = false;
                
            % In each trial, subsample spikes based on co-recorded msn
            % trial-wise spike count. If trial-wise spike count is < minTrialSpikes, or
            % the trial-wise spike count in all co-recorded MSNs is < minTrialSpikes,
            % skip that trial. If the number of FSI spikes in that trial <
            % the number of spikes in any co-recorded MSN for that trial,
            % don't subsamples and use all the FSI spikes
            for iS = 1:cfg_master.num_subsamples
                this_tw_ppc = nan(length(near_data.trial), length(od.fsi_res.near_spec{iM}.freqs));
                for iT = 1:length(near_data.trial)
                    valid_choices = od.msn_near_dist(iT,:);
                    valid_choices = valid_choices(valid_choices>=cfg_master.minTrialSpikes);
                    if od.fsi_res.near_spec{iM}.trialwise_spk_count(iT) < cfg_in.minTrialSpikes || isempty(valid_choices)
                        continue;
                    elseif od.fsi_res.near_spec{iM}.trialwise_spk_count(iT) <= max(od.msn_near_dist(iT,:)) % at least one co-recorded MSN has more spikes in this trial
                        % Don't subsample, use all spikes
                        sub_idx = trial_wise_spike{iT};
                        sub_sts = this_sts;
                        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                        this_tw_ppc(iT,:) = sub_ppc.ppc0;
                    else
                        % Choose a number to subsample based on co-recorded MSN spike counts for this trial
                        s_factor = 1/length(valid_choices);
                        choice = floor(rand()/s_factor)+1;                 
                        sub_idx = randsample(trial_wise_spike{iT}, valid_choices(choice));
                        sub_idx = sort(sub_idx);
                        sub_sts = this_sts;
                        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                        this_tw_ppc(iT,:) = sub_ppc.ppc0;
                    end
                end
                od.fsi_res.near_spec{iM}.trialwise_ppc(iS,:,:) = this_tw_ppc;
            end
        end 
    end
    if cfg_master.write_output  
         [~, fp, ~] = fileparts(pwd);
         pushdir(cfg_master.output_dir);
         fn_out = cat(2, fp, '_ft_spec.mat');
         save(fn_out,'od'); % should add option to save in specified output dir
         popdir;
    end
end

%% Other functions

%%
function S = LoadSpikesTarget(cfg_in)
    if ~isfield(cfg_in, 'Target') % no target specified, load them all
        S = LoadSpikes([]);
        return;
    end
    LoadExpKeys;
    % target specified, need to do some work
    % first see if this session has more than one target
    nTargets = length(ExpKeys.Target);
    if ~iscell(ExpKeys.Target) || (iscell(ExpKeys.Target) && length(ExpKeys.Target) == 1) % one target
        target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
        if isempty(target_idx)
            S = ts;
        else
            S = LoadSpikes([]);
        end
    else % multiple targets, assume TetrodeTargets exists
        please = []; please.getTTnumbers = 1;
        S = LoadSpikes(please);
        target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
        tt_num_keep = find(ExpKeys.TetrodeTargets == target_idx);
        keep = ismember(S.usr.tt_num, tt_num_keep);
        S = SelectTS([], S, keep);
    end
end