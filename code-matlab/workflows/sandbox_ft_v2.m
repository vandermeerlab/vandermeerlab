%% script to generate STA spectra and average STS on a trial-by trial basis as well as for binned trials
% trials are binned on the basis of mean firing rate in the trial such that
% the number of spikes in each of the bins are as close as possible
% The spectra for the FSIs are calculated after multiple rounds of
% subsampling
%% setup
clear;
cd('D:\ADRLabData');
% cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData');
please = [];
please.rats = {'R117', 'R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.output_dir = 'D:\RandomVstrAnalysis\temp';
% cfg_in.output_dir = '/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp';
cfg_in.incl_types = [1, 2];
cfg_in.nMinSpikes1 = 400; % For on track
cfg_in.nMinSpikes2 = 200; % For near and away
cfg_in.nMinSpikes3 = 100; % For lfr, hfr, p1 and p2


%%
% Top level loop which calls the main function for all the sessions
diary on;
for iS = 1:length(cfg_in.fd) % for each session...
    cfg_in.iS = iS;
    pushdir(cfg_in.fd{iS});
    generateSTS(cfg_in); % do the business
    popdir;
    
end % of sessions
diary off;

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
        % Restrict data to only on-track data
        on_track_data = ft_redefinetrial(cfg_onTrack, this_data);
        this_flag = false;
        % Sanity check to ensure that redefine trial has the same number of
        % spikes as restrict()
        if (sum(on_track_data.trial{1}(2,:)) ~= length(sd.S.t{iC}))
           this_flag = true;
           warning('ft_redefinetrial has %d spikes but restrict shows %d spikes on Track', ...
               sum(on_track_data.trial{1}(2,:)), length(sd.S.t{iC}))
        end
        od.msn_res.onTrack_spec{iM}.spk_count = sum(on_track_data.trial{1}(2,:));
        this_sta = ft_spiketriggeredaverage(cfg_ft, on_track_data);
        od.msn_res.onTrack_spec{iM}.sta_time = this_sta.time;
        od.msn_res.onTrack_spec{iM}.sta_vals = this_sta.avg(:,:)';
        od.msn_res.onTrack_spec{iM}.flag_unequalSpikes = this_flag;
        
        % Calculate and save STS
        cfg_sts.method = 'mtmconvol';
        cfg_sts.foi = 1:1:100;
        cfg_sts.t_ftimwin = 5./cfg_sts.foi;
        cfg_sts.taper = 'hanning';
        cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
        cfg_sts.channel = on_track_data.label{1};
        this_sts = ft_spiketriggeredspectrum(cfg_sts, on_track_data);
        this_flag = false;
        % Display warning to show that there were Nans in this calculation
        if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
            this_flag = true;
            warning('Cell %s has nans in its STS',sd.S.label{iC});
        end
        od.msn_res.onTrack_spec{iM}.freqs = this_sts.freq;
        od.msn_res.onTrack_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
        od.msn_res.onTrack_spec{iM}.flag_nansts = this_flag;
        
        % Calculate and save PPC
        cfg_ppc               = [];
        cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
        cfg_ppc.spikechannel  = this_sts.label;
        cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
        cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
        cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
        this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
        this_flag = false;
        % Display warning to show that there were Nans in this calculation
        if ~isempty(find(isnan(this_ppc.ppc0),1))
            this_flag = true;
            warning('Cell %s has nans in its ppc',sd.S.label{iC});
        end
        od.msn_res.onTrack_spec{iM}.ppc = this_ppc.ppc0';
        od.msn_res.onTrack_spec{iM}.flag_nanppc = this_flag;
        
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
        w_start1 = rt1 - 5;
        w_end1 = rt1 + 5;
        w_start2 = rt2 - 5;
        w_end2 =  rt2 + 5;
        % Last trial time shouldn't exceed Experiment end time
        w_end2(end) = min(w_end2(end), ExpKeys.TimeOffTrack);
        w_start = sort([w_start1; w_start2]);
        w_end = sort([w_end1; w_end2]);
        rt_iv = iv(w_start, w_end);
        % TODO: Catch warnings from Merge IV and set a flag
        rt_iv = MergeIV([], rt_iv);
        
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
        for iT = 1:length(near_data.trial)
           spk_count2 = spk_count2 + sum(near_data.trial{iT}(2,:)); 
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
            % Calculate and save STA
            this_sta = ft_spiketriggeredaverage(cfg_ft, near_data);
            od.msn_res.near_spec{iM}.sta_time = this_sta.time;
            od.msn_res.near_spec{iM}.sta_vals = this_sta.avg(:,:)';
            od.msn_res.near_spec{iM}.flag_unequalSpikes = this_flag;

            % Calculate and save STS
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

            % Calculate and save PPC
            cfg_ppc               = [];
            cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
            cfg_ppc.spikechannel  = this_sts.label;
            cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
            cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
            cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
            this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
            this_flag = false;
            % Display warning to show that there were Nans in this calculation
            if ~isempty(find(isnan(this_ppc.ppc0),1))
                this_flag = true;
                warning('Cell %s has nans in its ppc',sd.S.label{iC});
            end
            od.msn_res.near_spec{iM}.ppc = this_ppc.ppc0';
            od.msn_res.near_spec{iM}.flag_nanppc = this_flag;
            
            % Block of code to split trials into HFR and LFR
            tcount = length(near_data.trial);
            mfr = zeros(tcount, 1);
            all_tspikes = cell(1,tcount);
            for iT = 1:tcount
                all_tspikes{iT} = sum(near_data.trial{iT}(2,:));
                mfr(iT) = all_tspikes{iT}/near_data.time{iT}(end); 
            end
            % Get rid of trials with no spikes and bin the rest
            nz_trials = find(mfr ~= 0);
            nz_mfr = mfr(nz_trials);
            nz_tcount = length(nz_trials);
            nz_tspikes = cell(nz_tcount,1);
            spk_tcount = zeros(nz_tcount,1);
            for iT = 1:nz_tcount
                nz_tspikes{iT} = all_tspikes{nz_trials(iT)};
                spk_tcount(iT) = nz_tspikes{iT};
            end
            % Find out the firing rate threhsold to split the trials such that
            % spikes are more or less equally divided
            ufr = unique(nz_mfr);
            dif_min = sum(cell2mat(nz_tspikes));
            fr_thresh = 0;
            for iF = 1:length(ufr)
                cur_thresh = ufr(iF);
                l_spikes = sum(spk_tcount(nz_mfr <= cur_thresh));
                h_spikes = sum(spk_tcount(nz_mfr > cur_thresh));
                cur_dif = abs(h_spikes - l_spikes);
                if dif_min > cur_dif
                    dif_min = cur_dif;
                    fr_thresh = cur_thresh;
                end
            end            
            hfr_trials = mfr > fr_thresh;
            lfr_trials = ~hfr_trials;
            
            
            % Divide trials into partitions with equal spikes but not on
            % the basis of fr_threshold
            % make 2 groups and put the largest remaining element to the partition
            % with the smallest sum. Once a sum reaches n/2, put the rest
            % of the trials in the other group
            % Potential Problem: The remaining lower spike might be the
            % lfr trials, and thus biasing the partition.
            total_spks = sum(spk_tcount);
            spk_thresh = total_spks/2;
            [sorted_tcount, sorted_idx] = sort(cell2mat(all_tspikes), 'descend');
            trial_part = zeros(size(sorted_tcount));
            g1_sum = 0;
            g2_sum = 0;
            for iT = 1:length(sorted_tcount)
                if g1_sum <= g2_sum
                    trial_part(iT) = 1;
                else
                    trial_part(iT) = 2;
                end
                g1_sum = sum(sorted_tcount(trial_part == 1));
                g2_sum = sum(sorted_tcount(trial_part == 2));
                if g1_sum >= spk_thresh
                    trial_part(iT+1:end) = 2;
                    break;
                end
                if g2_sum >= spk_thresh
                    trial_part(iT+1:end) = 1;
                    break;
                end
            end
            p1_trials = false(size(sorted_tcount));
            p1_trials(sorted_idx(trial_part==1)) = true;
            p2_trials = ~p1_trials;
           
            od.msn_res.near_spec{iM}.mfr = mfr;
            od.msn_res.near_spec{iM}.fr_thresh = fr_thresh;
            od.msn_res.near_spec{iM}.trial_spk_count = cell2mat(all_tspikes);
            od.msn_res.near_spec{iM}.p1_trials = p1_trials;
            od.msn_res.near_spec{iM}.p2_trials = p2_trials;

            % Calculate and Save all spec results for near_hfr_data
            cfg_near_hfr_trials.trl = cfg_near_trials.trl(hfr_trials,:);
            near_hfr_data = ft_redefinetrial(cfg_near_hfr_trials, this_data);
            
            % Skip if no spks in near_hfr
            spk_count = 0;
            for iT = 1:length(near_hfr_data.trial)
                spk_count = spk_count + sum(near_hfr_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes2
                od.msn_res.near_hfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.near_hfr_spec{iM}.spk_count = spk_count;
                od.msn_res.near_hfr_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, near_hfr_data);
                od.msn_res.near_hfr_spec{iM}.sta_time = this_sta.time;
                od.msn_res.near_hfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = near_hfr_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, near_hfr_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its near STS',sd.S.label{iC});
                end
                od.msn_res.near_hfr_spec{iM}.flag_nansts = this_flag;
                od.msn_res.near_hfr_spec{iM}.freqs = this_sts.freq;
                od.msn_res.near_hfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));

                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.near_hfr_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.near_hfr_spec{iM}.flag_nanppc = this_flag;
            end      
                       
            % Calculate and Save all spec results for near_lfr_data
            cfg_near_lfr_trials.trl = cfg_near_trials.trl(lfr_trials,:);
            near_lfr_data = ft_redefinetrial(cfg_near_lfr_trials, this_data);
            
            % Skip if no spks in near_lfr
            spk_count = 0;
            for iT = 1:length(near_lfr_data.trial)
                spk_count = spk_count + sum(near_lfr_data.trial{iT}(2,:)); 
            end
            if spk_count <  cfg_master.nMinSpikes3
                od.msn_res.near_lfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.near_lfr_spec{iM}.spk_count = spk_count;
                od.msn_res.near_lfr_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, near_lfr_data);
                od.msn_res.near_lfr_spec{iM}.sta_time = this_sta.time;
                od.msn_res.near_lfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = near_lfr_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, near_lfr_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its STS',sd.S.label{iC});
                end
                od.msn_res.near_lfr_spec{iM}.freqs = this_sts.freq;
                od.msn_res.near_lfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                od.msn_res.near_lfr_spec{iM}.flag_nansts = this_flag;

                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.near_lfr_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.near_lfr_spec{iM}.flag_nanppc = this_flag;   
            end
            
            % Calculate and Save all spec results for near_p1_data
            cfg_near_p1_trials.trl = cfg_near_trials.trl(p1_trials,:);
            near_p1_data = ft_redefinetrial(cfg_near_p1_trials, this_data);
            
            % Skip if no spks in near_p1
            spk_count = 0;
            for iT = 1:length(near_p1_data.trial)
                spk_count = spk_count + sum(near_p1_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes3
                od.msn_res.near_p1_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.near_p1_spec{iM}.spk_count = spk_count;
                od.msn_res.near_p1_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, near_p1_data);
                od.msn_res.near_p1_spec{iM}.sta_time = this_sta.time;
                od.msn_res.near_p1_spec{iM}.sta_vals = this_sta.avg(:,:)';

                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = near_p1_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, near_p1_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its near STS',sd.S.label{iC});
                end
                od.msn_res.near_p1_spec{iM}.flag_nansts = this_flag;
                od.msn_res.near_p1_spec{iM}.freqs = this_sts.freq;
                od.msn_res.near_p1_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));

                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.near_p1_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.near_p1_spec{iM}.flag_nanppc = this_flag;
            end
            
            % Calculate and Save all spec results for near_p2_data
            cfg_near_p2_trials.trl = cfg_near_trials.trl(p2_trials,:);
            near_p2_data = ft_redefinetrial(cfg_near_p2_trials, this_data);
            
            % Skip if no spks in near_p2
            spk_count = 0;
            for iT = 1:length(near_p2_data.trial)
                spk_count = spk_count + sum(near_p2_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes3
                od.msn_res.near_p2_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.near_p2_spec{iM}.spk_count = spk_count;
                od.msn_res.near_p2_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, near_p2_data);
                od.msn_res.near_p2_spec{iM}.sta_time = this_sta.time;
                od.msn_res.near_p2_spec{iM}.sta_vals = this_sta.avg(:,:)';

                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = near_p2_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, near_p2_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its near STS',sd.S.label{iC});
                end
                od.msn_res.near_p2_spec{iM}.flag_nansts = this_flag;
                od.msn_res.near_p2_spec{iM}.freqs = this_sts.freq;
                od.msn_res.near_p2_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));

                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.near_p2_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.near_p2_spec{iM}.flag_nanppc = this_flag;
            end     
        end 
        
        % For away_reward_trials
        w_start = [ExpKeys.TimeOnTrack;rt2(1:end-1)+5];
        w_end = (rt1-5);
        % Get rid of extra long-trials (greater than mean + 1*SD)
        trial_length = w_end - w_start;
        mtl = mean(trial_length); stl = std(trial_length);
        valid_trials = (trial_length <= mtl + stl); 
        w_start = w_start(valid_trials);
        w_end= w_end(valid_trials);
        rt_iv = iv(w_start, w_end);
        % TODO: Catch warnings from Merge IV and set a flag
        rt_iv = MergeIV([], rt_iv);
        
        % Break down data into away trials
        temp_tvec = ft_csc.time{1} + double(ft_csc.hdr.FirstTimeStamp)/1e6;     
        temp_start = nearest_idx3(rt_iv.tstart, temp_tvec);
        temp_end = nearest_idx3(rt_iv.tend, temp_tvec);
        cfg_away_trials.trl = [temp_start, temp_end, zeros(size(temp_start))];
        away_data = ft_redefinetrial(cfg_away_trials, this_data);
        % Sanity check to ensure that redefine trial has the same number of
        % spikes as restrict()
        this_flag = false;
        spk_count1 = length(restrict(sd.S, rt_iv).t{iC});
        spk_count2 = 0;
        for iT = 1:length(away_data.trial)
           spk_count2 = spk_count2 + sum(away_data.trial{iT}(2,:)); 
        end
        % Skip if no spikes present!
        if spk_count2 < cfg_master.nMinSpikes2
            od.msn_res.away_spec{iM}.flag_tooFewSpikes = true;
        else
            od.msn_res.away_spec{iM}.spk_count = spk_count2;
            od.msn_res.away_spec{iM}.flag_tooFewSpikes = false;
            if spk_count1 ~= spk_count2
                this_flag = true;
                warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in away-Reward trials', ...
                    spk_count1, spk_count2)

            end
            % Calculate and save STA
            this_sta = ft_spiketriggeredaverage(cfg_ft, away_data);
            od.msn_res.away_spec{iM}.sta_time = this_sta.time;
            od.msn_res.away_spec{iM}.sta_vals = this_sta.avg(:,:)';
            od.msn_res.away_spec{iM}.flag_unequalSpikes = this_flag;

            % Calculate and save STS
            cfg_sts.method = 'mtmconvol';
            cfg_sts.foi = 1:1:100;
            cfg_sts.t_ftimwin = 5./cfg_sts.foi;
            cfg_sts.taper = 'hanning';
            cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
            cfg_sts.channel = away_data.label{1};
            this_sts = ft_spiketriggeredspectrum(cfg_sts, away_data);
            this_flag = false;
            % Display warning to show that there were Nans in this calculation
            if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                this_flag = true;
                warning('Cell %s has nans in its STS',sd.S.label{iC});
            end
            od.msn_res.away_spec{iM}.freqs = this_sts.freq;
            od.msn_res.away_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
            od.msn_res.away_spec{iM}.flag_nansts = this_flag;

            % Calculate and save PPC
            cfg_ppc               = [];
            cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
            cfg_ppc.spikechannel  = this_sts.label;
            cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
            cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
            cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
            this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
            this_flag = false;
            % Display warning to show that there were Nans in this calculation
            if ~isempty(find(isnan(this_ppc.ppc0),1))
                this_flag = true;
                warning('Cell %s has nans in its ppc',sd.S.label{iC});
            end
            od.msn_res.away_spec{iM}.ppc = this_ppc.ppc0';
            od.msn_res.away_spec{iM}.flag_nanppc = this_flag;
            
            % Block of code to split trials into HFR and LFR
            tcount = length(away_data.trial);
            mfr = zeros(tcount, 1);
            all_tspikes = cell(1,tcount);
            for iT = 1:tcount
                all_tspikes{iT} = sum(away_data.trial{iT}(2,:));
                mfr(iT) = all_tspikes{iT}/away_data.time{iT}(end); 
            end
            % Get rid of trials with no spikes and bin the rest
            nz_trials = find(mfr ~= 0);
            nz_mfr = mfr(nz_trials);
            nz_tcount = length(nz_trials);
            nz_tspikes = cell(nz_tcount,1);
            spk_tcount = zeros(nz_tcount,1);
            for iT = 1:nz_tcount
                nz_tspikes{iT} = all_tspikes{nz_trials(iT)};
                spk_tcount(iT) = nz_tspikes{iT};
            end
            % Find out the firing rate threhsold to split the trials such that
            % spikes are more or less equally divided
            ufr = unique(nz_mfr);
            dif_min = sum(cell2mat(nz_tspikes));
            fr_thresh = 0;
            for iF = 1:length(ufr)
                cur_thresh = ufr(iF);
                l_spikes = sum(spk_tcount(nz_mfr <= cur_thresh));
                h_spikes = sum(spk_tcount(nz_mfr > cur_thresh));
                cur_dif = abs(h_spikes - l_spikes);
                if dif_min > cur_dif
                    dif_min = cur_dif;
                    fr_thresh = cur_thresh;
                end
            end            
            hfr_trials = mfr > fr_thresh;
            lfr_trials = ~hfr_trials;
            
            % Divide trials into partitions with equal spikes but not on
            % the basis of fr_threshold
            % make 2 groups and put the largest remaining element to the partition
            % with the smallest sum. Once a sum reaches n/2, put the rest
            % of the trials in the other group
            % Potential Problem: The remaining lower spike might be the
            % lfr trials, and thus biasing the partition.
            total_spks = sum(spk_tcount);
            spk_thresh = total_spks/2;
            [sorted_tcount, sorted_idx] = sort(cell2mat(all_tspikes), 'descend');
            trial_part = zeros(size(sorted_tcount));
            g1_sum = 0;
            g2_sum = 0;
            for iT = 1:length(sorted_tcount)
                if g1_sum <= g2_sum
                    trial_part(iT) = 1;
                else
                    trial_part(iT) = 2;
                end
                g1_sum = sum(sorted_tcount(trial_part == 1));
                g2_sum = sum(sorted_tcount(trial_part == 2));
                if g1_sum >= spk_thresh
                    trial_part(iT+1:end) = 2;
                    break;
                end
                if g2_sum >= spk_thresh
                    trial_part(iT+1:end) = 1;
                    break;
                end
            end
            p1_trials = false(size(sorted_tcount));
            p1_trials(sorted_idx(trial_part==1)) = true;
            p2_trials = ~p1_trials;
           
            od.msn_res.away_spec{iM}.mfr = mfr;
            od.msn_res.away_spec{iM}.fr_thresh = fr_thresh;
            od.msn_res.away_spec{iM}.trial_spk_count = cell2mat(all_tspikes);
            od.msn_res.away_spec{iM}.p1_trials = p1_trials;
            od.msn_res.away_spec{iM}.p2_trials = p2_trials;

            % Calculate and Save all spec results for away_hfr_data
            cfg_away_hfr_trials.trl = cfg_away_trials.trl(hfr_trials,:);
            away_hfr_data = ft_redefinetrial(cfg_away_hfr_trials, this_data);
            
            % Skip if no spks in away_hfr
            spk_count = 0;
            for iT = 1:length(away_hfr_data.trial)
                spk_count = spk_count + sum(away_hfr_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes3
                od.msn_res.away_hfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.away_hfr_spec{iM}.spk_count = spk_count;
                od.msn_res.away_hfr_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, away_hfr_data);
                od.msn_res.away_hfr_spec{iM}.sta_time = this_sta.time;
                od.msn_res.away_hfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = away_hfr_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, away_hfr_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its away STS',sd.S.label{iC});
                end
                od.msn_res.away_hfr_spec{iM}.flag_nansts = this_flag;
                od.msn_res.away_hfr_spec{iM}.freqs = this_sts.freq;
                od.msn_res.away_hfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));

                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.away_hfr_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.away_hfr_spec{iM}.flag_nanppc = this_flag;
            end      
                       
            % Calculate and Save all spec results for away_lfr_data
            cfg_away_lfr_trials.trl = cfg_away_trials.trl(lfr_trials,:);
            away_lfr_data = ft_redefinetrial(cfg_away_lfr_trials, this_data);
            
            % Skip if no spks in away_lfr
            spk_count = 0;
            for iT = 1:length(away_lfr_data.trial)
                spk_count = spk_count + sum(away_lfr_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes3
                od.msn_res.away_lfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.away_lfr_spec{iM}.spk_count = spk_count;
                od.msn_res.away_lfr_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, away_lfr_data);
                od.msn_res.away_lfr_spec{iM}.sta_time = this_sta.time;
                od.msn_res.away_lfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = away_lfr_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, away_lfr_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its STS',sd.S.label{iC});
                end
                od.msn_res.away_lfr_spec{iM}.freqs = this_sts.freq;
                od.msn_res.away_lfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                od.msn_res.away_lfr_spec{iM}.flag_nansts = this_flag;

                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.away_lfr_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.away_lfr_spec{iM}.flag_nanppc = this_flag;   
            end
            
            % Calculate and Save all spec results for away_p1_data
            cfg_away_p1_trials.trl = cfg_away_trials.trl(p1_trials,:);
            away_p1_data = ft_redefinetrial(cfg_away_p1_trials, this_data);
            
            % Skip if no spks in away_p1
            spk_count = 0;
            for iT = 1:length(away_p1_data.trial)
                spk_count = spk_count + sum(away_p1_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes3
                od.msn_res.away_p1_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.away_p1_spec{iM}.spk_count = spk_count;
                od.msn_res.away_p1_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, away_p1_data);
                od.msn_res.away_p1_spec{iM}.sta_time = this_sta.time;
                od.msn_res.away_p1_spec{iM}.sta_vals = this_sta.avg(:,:)';

                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = away_p1_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, away_p1_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its away STS',sd.S.label{iC});
                end
                od.msn_res.away_p1_spec{iM}.flag_nansts = this_flag;
                od.msn_res.away_p1_spec{iM}.freqs = this_sts.freq;
                od.msn_res.away_p1_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));

                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.away_p1_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.away_p1_spec{iM}.flag_nanppc = this_flag;
            end
            
            % Calculate and Save all spec results for away_p2_data
            cfg_away_p2_trials.trl = cfg_away_trials.trl(p2_trials,:);
            away_p2_data = ft_redefinetrial(cfg_away_p2_trials, this_data);
            
            % Skip if no spks in away_p2
            spk_count = 0;
            for iT = 1:length(away_p2_data.trial)
                spk_count = spk_count + sum(away_p2_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes3
                od.msn_res.away_p2_spec{iM}.flag_tooFewSpikes = true;
            else
                od.msn_res.away_p2_spec{iM}.spk_count = spk_count;
                od.msn_res.away_p2_spec{iM}.flag_tooFewSpikes = false;
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, away_p2_data);
                od.msn_res.away_p2_spec{iM}.sta_time = this_sta.time;
                od.msn_res.away_p2_spec{iM}.sta_vals = this_sta.avg(:,:)';
                % Calculate and save STS
                cfg_sts.method = 'mtmconvol';
                cfg_sts.foi = 1:1:100;
                cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                cfg_sts.taper = 'hanning';
                cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                cfg_sts.channel = away_p2_data.label{1};
                this_sts = ft_spiketriggeredspectrum(cfg_sts, away_p2_data);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_sts.fourierspctrm{1}),1))
                    this_flag = true;
                    warning('Cell %s has nans in its away STS',sd.S.label{iC});
                end
                od.msn_res.away_p2_spec{iM}.flag_nansts = this_flag;
                od.msn_res.away_p2_spec{iM}.freqs = this_sts.freq;
                od.msn_res.away_p2_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                % Calculate and save PPC
                cfg_ppc               = [];
                cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                cfg_ppc.spikechannel  = this_sts.label;
                cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                this_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_sts);
                this_flag = false;
                % Display warning to show that there were Nans in this calculation
                if ~isempty(find(isnan(this_ppc.ppc0),1))
                    this_flag = true;
                    warning('Cell %s has nans in its ppc',sd.S.label{iC});
                end
                od.msn_res.away_p2_spec{iM}.ppc = this_ppc.ppc0';
                od.msn_res.away_p2_spec{iM}.flag_nanppc = this_flag;
            end     
        end
    end
    
    % Get the msn distributions
    msn_onTrack_dist = [];
    msn_near_dist = [];
    msn_near_lfr_dist = [];
    msn_near_hfr_dist = [];
    msn_near_p1_dist = [];
    msn_near_p2_dist = [];
    msn_away_dist = [];
    msn_away_lfr_dist = [];
    msn_away_hfr_dist = [];
    msn_away_p1_dist = [];
    msn_away_p2_dist = [];

    for iM = 1:length(od.msn_res.onTrack_spec)
        msn_onTrack_dist = [msn_onTrack_dist od.msn_res.onTrack_spec{iM}.spk_count];
        if ~od.msn_res.near_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_spec{iM}.flag_nansts & ~od.msn_res.near_spec{iM}.flag_nanppc
            msn_near_dist = [msn_near_dist od.msn_res.near_spec{iM}.spk_count];
            if ~od.msn_res.near_lfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_lfr_spec{iM}.flag_nansts & ~od.msn_res.near_lfr_spec{iM}.flag_nanppc
                msn_near_lfr_dist = [msn_near_lfr_dist od.msn_res.near_lfr_spec{iM}.spk_count];
            end
            if ~od.msn_res.near_hfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_hfr_spec{iM}.flag_nansts & ~od.msn_res.near_hfr_spec{iM}.flag_nanppc
                msn_near_hfr_dist = [msn_near_hfr_dist od.msn_res.near_hfr_spec{iM}.spk_count];
            end
            if ~od.msn_res.near_p1_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_p1_spec{iM}.flag_nansts & ~od.msn_res.near_p1_spec{iM}.flag_nanppc
                msn_near_p1_dist = [msn_near_p1_dist od.msn_res.near_p1_spec{iM}.spk_count];
            end
            if ~od.msn_res.near_p2_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_p2_spec{iM}.flag_nansts & ~od.msn_res.near_p2_spec{iM}.flag_nanppc
                msn_near_p2_dist = [msn_near_p2_dist od.msn_res.near_p2_spec{iM}.spk_count];
            end
        end
        if ~od.msn_res.away_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_spec{iM}.flag_nansts & ~od.msn_res.away_spec{iM}.flag_nanppc
            msn_away_dist = [msn_away_dist od.msn_res.away_spec{iM}.spk_count];
            if ~od.msn_res.away_lfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_lfr_spec{iM}.flag_nansts & ~od.msn_res.away_lfr_spec{iM}.flag_nanppc
                msn_away_lfr_dist = [msn_away_lfr_dist od.msn_res.away_lfr_spec{iM}.spk_count];
            end
            if ~od.msn_res.away_hfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_hfr_spec{iM}.flag_nansts & ~od.msn_res.away_hfr_spec{iM}.flag_nanppc
                msn_away_hfr_dist = [msn_away_hfr_dist od.msn_res.away_hfr_spec{iM}.spk_count];
            end
            if ~od.msn_res.away_p1_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_p1_spec{iM}.flag_nansts & ~od.msn_res.away_p1_spec{iM}.flag_nanppc
                msn_away_p1_dist = [msn_away_p1_dist od.msn_res.away_p1_spec{iM}.spk_count];
            end
            if ~od.msn_res.away_p2_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_p2_spec{iM}.flag_nansts & ~od.msn_res.away_p2_spec{iM}.flag_nanppc
                msn_away_p2_dist = [msn_away_p2_dist od.msn_res.away_p2_spec{iM}.spk_count];
            end
        end
    end
    
    all_fsi = find(od.S1.cell_type == 1);

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

