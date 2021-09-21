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
cfg_in.nMinSpikes2 = 400; % For near and away
cfg_in.nMinSpikes3 = 200; % For lfr, hfr, p1 and p2
cfg_in.nControlSplits = 2;
cfg_in.num_subsamples = 2;


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
            % Find out the firing rate threshold to split the trials such that
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
            
            
            % Divide trials into nControlSplits partitions with nearly equal spikes but not on
            % the basis of fr_threshold. 
            % Algo used: Use randomly generated splits and accept a split only if
            % 1) difference between splits is <= 2x the diffference of the
            % FR based split
            % 2) not the same as an experimental split
            % 3) not the same as an already found split
            % Wait for 100,000,000 iterations, if you haven't found
            % nControlSplits valid splits by then, take what you have!
            
            valid_splits  = false(cfg_master.nControlSplits,tcount);
            all_tspikes_mat = cell2mat(all_tspikes);
            last_valid_split = 0;
            for iRand = 1:100000000 % wait till 100 million iterations
                A = 1:tcount;
                ndiv = 2;
                this_idx = sort([1 randperm(length(A)-1, ndiv-1)+1 length(A)+1]);
                for k1 = 1:length(this_idx)-1
                    R{k1} = A(this_idx(k1):this_idx(k1+1)-1);
                end  
                this_perm = randperm(tcount);
                this_split = false(1,tcount);
                this_split(this_perm(R{1})) = true;
                this_split_dif = abs(sum(all_tspikes_mat(this_split)) - sum(all_tspikes_mat(~this_split)));
                % reject split of split dif is grater than threshold
                if this_split_dif > dif_min
                   continue; 
                end
                % reject if split exactly the same as hypothesis split
                if sum(this_split == hfr_trials) == tcount | sum(this_split == lfr_trials) == tcount
                   continue;
                end
                % reject if split the same as previously found valid split
                flag_repeat_split = false;
                for iCheck = 1:1:last_valid_split
                    if sum(this_split == valid_splits(iCheck)) == tcount | sum(~this_split == valid_splits(iCheck)) == tcount
                        flag_repeat_spilt = true;
                        break;
                    end
                end
                if flag_repeat_split
                    continue;
                end
                % If you have made it till here, you found a valid_split!
                last_valid_split = last_valid_split + 1;
                valid_splits(last_valid_split,:) = this_split;
                % If 100 valid splits are found, get out of this!
                if last_valid_split == cfg_master.nControlSplits
                    break;
                end
            end
              
            od.msn_res.near_spec{iM}.mfr = mfr;
            od.msn_res.near_spec{iM}.fr_thresh = fr_thresh;
            od.msn_res.near_spec{iM}.trial_spk_count = cell2mat(all_tspikes);
            od.msn_res.near_spec{iM}.valid_split_count = last_valid_split;
            od.msn_res.near_spec{iM}.valid_splits = valid_splits;
            od.msn_res.near_spec{iM}.randIters = iRand;

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
            
            % Do control split stuff only if lfr/mfr splits are perfect
            % i.e., enough spikes and no nan flags
            
            if ~od.msn_res.near_lfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.msn_res.near_hfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.msn_res.near_lfr_spec{iM}.flag_nansts && ...
                    ~od.msn_res.near_lfr_spec{iM}.flag_nanppc  && ...
                    ~od.msn_res.near_hfr_spec{iM}.flag_nansts && ...
                    ~od.msn_res.near_hfr_spec{iM}.flag_nanppc && ...
                    od.msn_res.near_spec{iM}.valid_split_count == cfg_master.nControlSplits
                
                od.msn_res.near_spec{iM}.flag_no_control_split = false;
                
                % Create temporary arrays to average later
                all_sta_vals = zeros([2, cfg_master.nControlSplits, size(od.msn_res.near_spec{iM}.sta_vals)]);
                all_freqs = zeros([cfg_master.nControlSplits, size(od.msn_res.near_spec{iM}.freqs)]);
                all_sts_vals = zeros([2, cfg_master.nControlSplits, length(od.msn_res.near_spec{iM}.sts_vals)]);
                all_ppc = zeros([2, cfg_master.nControlSplits, size(od.msn_res.near_spec{iM}.ppc)]);
                all_spk_count = zeros(2, cfg_master.nControlSplits);
         
                
                % For each split, fill in the tables
                % If any one of them is nan, exit loop, set
                % flag_no_control_split as true and move on to the next
                % cell    
                flag_nan_in_split = false;      
                for iSplit = 1:cfg_master.nControlSplits
  
                    % Calculate and Save all spec results for near_lfr_data
                    this_p1_trials = find(od.msn_res.near_spec{iM}.valid_splits(iSplit,:));
                    this_p1_cfg.trl = cfg_near_trials.trl(this_p1_trials,:);
                    this_p1_data = ft_redefinetrial(this_p1_cfg, this_data);
                    this_p2_trials = find(~od.msn_res.near_spec{iM}.valid_splits(iSplit,:));
                    this_p2_cfg.trl = cfg_near_trials.trl(this_p2_trials,:);
                    this_p2_data = ft_redefinetrial(this_p2_cfg, this_data);
                    
                    p1_spk_count = 0;
                    for iT = 1:length(this_p1_data.trial)
                        p1_spk_count = p1_spk_count + sum(this_p1_data.trial{iT}(2,:)); 
                    end
                    
                    p2_spk_count = 0;
                    for iT = 1:length(this_p2_data.trial)
                        p2_spk_count = p2_spk_count + sum(this_p2_data.trial{iT}(2,:)); 
                    end
                    
                    all_spk_count(1,iSplit) = p1_spk_count;
                    all_spk_count(2,iSplit) = p2_spk_count;

                    % Calculate and save STA
                    this_p1_sta = ft_spiketriggeredaverage(cfg_ft, this_p1_data);
                    this_p2_sta = ft_spiketriggeredaverage(cfg_ft, this_p2_data);
                    all_sta_vals(1, iSplit, :) = this_p1_sta.avg(:,:)';
                    all_sta_vals(2, iSplit, :) = this_p2_sta.avg(:,:)';
                    
                    % Calculate and save STS
                    cfg_sts.method = 'mtmconvol';
                    cfg_sts.foi = 1:1:100;
                    cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                    cfg_sts.taper = 'hanning';
                    cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                    cfg_sts.channel = near_lfr_data.label{1};
                    this_p1_sts = ft_spiketriggeredspectrum(cfg_sts, this_p1_data);
                    this_p2_sts = ft_spiketriggeredspectrum(cfg_sts, this_p2_data);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_sts.fourierspctrm{1}),1)) || ...
                            ~isempty(find(isnan(this_p2_sts.fourierspctrm{1}),1))
                        flag_nan_in_split = true;
                        break;
                    end
                    all_freqs(iSplit,:) = this_p1_sts.freq;
                    all_sts_vals(1,iSplit,:) = nanmean(sq(abs(this_p1_sts.fourierspctrm{1})));
                    all_sts_vals(2,iSplit,:) = nanmean(sq(abs(this_p2_sts.fourierspctrm{1})));
                    
                    % Calculate and save PPC
                    cfg_ppc               = [];
                    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg_ppc.spikechannel  = this_sts.label;
                    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                    cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                    cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                    this_p1_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p1_sts);
                    this_p2_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p2_sts);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_ppc.ppc0),1)) || ...
                         ~isempty(find(isnan(this_p2_ppc.ppc0),1))   
                        flag_nan_in_split = true;
                        break;
                    end
                    all_ppc(1,iSplit,:) = this_p1_ppc.ppc0';
                    all_ppc(2,iSplit,:) = this_p2_ppc.ppc0'; 
                end
                if flag_nan_in_split
                   od.msn_res.near_spec{iM}.flag_no_control_split = true;
                   break;
                end
                
                % Save averages and std
                od.msn_res.near_p1_spec{iM}.sta_time = od.msn_res.near_spec{iM}.sta_time;
                od.msn_res.near_p2_spec{iM}.sta_time = od.msn_res.near_spec{iM}.sta_time;
                od.msn_res.near_p1_spec{iM}.freqs = squeeze(mean(all_freqs(1,:,:),2))';
                od.msn_res.near_p2_spec{iM}.freqs = squeeze(mean(all_freqs(2,:,:),2))';
                od.msn_res.near_p1_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(1,:,:),2))';
                od.msn_res.near_p2_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(2,:,:),2))';
                od.msn_res.near_p1_spec{iM}.sd_sta = squeeze(std(all_sta_vals(1,:,:)))';
                od.msn_res.near_p2_spec{iM}.sd_sta = squeeze(std(all_sta_vals(2,:,:)))';
                od.msn_res.near_p1_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(1,:,:),2))';
                od.msn_res.near_p2_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(2,:,:),2))';
                od.msn_res.near_p1_spec{iM}.sd_sts = squeeze(std(all_sts_vals(1,:,:)))';
                od.msn_res.near_p2_spec{iM}.sd_sts = squeeze(std(all_sts_vals(2,:,:)))';
                od.msn_res.near_p1_spec{iM}.mean_ppc = squeeze(mean(all_ppc(1,:,:),2))';
                od.msn_res.near_p2_spec{iM}.mean_ppc = squeeze(mean(all_ppc(2,:,:),2))';
                od.msn_res.near_p1_spec{iM}.sd_ppc = squeeze(std(all_ppc(1,:,:)))';
                od.msn_res.near_p2_spec{iM}.sd_ppc = squeeze(std(all_ppc(2,:,:)))';
                od.msn_res.near_p1_spec{iM}.spk_count = round(mean(all_spk_count(1,:)));
                od.msn_res.near_p2_spec{iM}.spk_count = round(mean(all_spk_count(2,:)));
                       
            else
                od.msn_res.near_spec{iM}.flag_no_control_split = true;
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
            % Find out the firing rate threshold to split the trials such that
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
            
            % Divide trials into nControlSplits partitions with awayly equal spikes but not on
            % the basis of fr_threshold. 
            % Algo used: Use randomly generated splits and accept a split only if
            % 1) difference between splits is <= 2x the diffference of the
            % FR based split
            % 2) not the same as an experimental split
            % 3) not the same as an already found split
            % Wait for 100,000,000 iterations, if you haven't found
            % nControlSplits valid splits by then, take what you have!
            
            valid_splits  = false(cfg_master.nControlSplits,tcount);
            all_tspikes_mat = cell2mat(all_tspikes);
            last_valid_split = 0;
            for iRand = 1:100000000 % wait till 100 million iterations
                A = 1:tcount;
                ndiv = 2;
                this_idx = sort([1 randperm(length(A)-1, ndiv-1)+1 length(A)+1]);
                for k1 = 1:length(this_idx)-1
                    R{k1} = A(this_idx(k1):this_idx(k1+1)-1);
                end  
                this_perm = randperm(tcount);
                this_split = false(1,tcount);
                this_split(this_perm(R{1})) = true;
                this_split_dif = abs(sum(all_tspikes_mat(this_split)) - sum(all_tspikes_mat(~this_split)));
                % reject split of split dif is grater than threshold
                if this_split_dif > dif_min
                   continue; 
                end
                % reject if split exactly the same as hypothesis split
                if sum(this_split == hfr_trials) == tcount | sum(this_split == lfr_trials) == tcount
                   continue;
                end
                % reject if split the same as previously found valid split
                flag_repeat_split = false;
                for iCheck = 1:1:last_valid_split
                    if sum(this_split == valid_splits(iCheck)) == tcount | sum(~this_split == valid_splits(iCheck)) == tcount
                        flag_repeat_spilt = true;
                        break;
                    end
                end
                if flag_repeat_split
                    continue;
                end
                % If you have made it till here, you found a valid_split!
                last_valid_split = last_valid_split + 1;
                valid_splits(last_valid_split,:) = this_split;
                % If 100 valid splits are found, get out of this!
                if last_valid_split == cfg_master.nControlSplits
                    break;
                end
            end
              
            od.msn_res.away_spec{iM}.mfr = mfr;
            od.msn_res.away_spec{iM}.fr_thresh = fr_thresh;
            od.msn_res.away_spec{iM}.trial_spk_count = cell2mat(all_tspikes);
            od.msn_res.away_spec{iM}.valid_split_count = last_valid_split;
            od.msn_res.away_spec{iM}.valid_splits = valid_splits;
            od.msn_res.away_spec{iM}.randIters = iRand;

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
            
                        % Do control split stuff only if lfr/mfr splits are perfect
            % i.e., enough spikes and no nan flags
            
            if ~od.msn_res.away_lfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.msn_res.away_hfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.msn_res.away_lfr_spec{iM}.flag_nansts && ...
                    ~od.msn_res.away_lfr_spec{iM}.flag_nanppc  && ...
                    ~od.msn_res.away_hfr_spec{iM}.flag_nansts && ...
                    ~od.msn_res.away_hfr_spec{iM}.flag_nanppc && ...
                    od.msn_res.away_spec{iM}.valid_split_count == cfg_master.nControlSplits
                
                od.msn_res.away_spec{iM}.flag_no_control_split = false;
                
                % Create temporary arrays to average later
                all_sta_vals = zeros([2, cfg_master.nControlSplits, size(od.msn_res.away_spec{iM}.sta_vals)]);
                all_freqs = zeros([cfg_master.nControlSplits, size(od.msn_res.away_spec{iM}.freqs)]);
                all_sts_vals = zeros([2, cfg_master.nControlSplits, length(od.msn_res.away_spec{iM}.sts_vals)]);
                all_ppc = zeros([2, cfg_master.nControlSplits, size(od.msn_res.away_spec{iM}.ppc)]);
                all_spk_count = zeros(2, cfg_master.nControlSplits);
         
                % For each split, fill in the tables
                % If any one of them is nan, exit loop, set
                % flag_no_control_split as true and move on to the next
                % cell    
                flag_nan_in_split = false;      
                for iSplit = 1:cfg_master.nControlSplits
                    % Calculate and Save all spec results for away_lfr_data
                    this_p1_trials = find(od.msn_res.away_spec{iM}.valid_splits(iSplit,:));
                    this_p1_cfg.trl = cfg_away_trials.trl(this_p1_trials,:);
                    this_p1_data = ft_redefinetrial(this_p1_cfg, this_data);
                    this_p2_trials = find(~od.msn_res.away_spec{iM}.valid_splits(iSplit,:));
                    this_p2_cfg.trl = cfg_away_trials.trl(this_p2_trials,:);
                    this_p2_data = ft_redefinetrial(this_p2_cfg, this_data);
                    
                    p1_spk_count = 0;
                    for iT = 1:length(this_p1_data.trial)
                        p1_spk_count = p1_spk_count + sum(this_p1_data.trial{iT}(2,:)); 
                    end
                    
                    p2_spk_count = 0;
                    for iT = 1:length(this_p2_data.trial)
                        p2_spk_count = p2_spk_count + sum(this_p2_data.trial{iT}(2,:)); 
                    end
                    
                    all_spk_count(1,iSplit) = p1_spk_count;
                    all_spk_count(2,iSplit) = p2_spk_count;
                    
                    % Calculate and save STA
                    this_p1_sta = ft_spiketriggeredaverage(cfg_ft, this_p1_data);
                    this_p2_sta = ft_spiketriggeredaverage(cfg_ft, this_p2_data);
                    all_sta_vals(1, iSplit, :) = this_p1_sta.avg(:,:)';
                    all_sta_vals(2, iSplit, :) = this_p2_sta.avg(:,:)';
                    
                    % Calculate and save STS
                    cfg_sts.method = 'mtmconvol';
                    cfg_sts.foi = 1:1:100;
                    cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                    cfg_sts.taper = 'hanning';
                    cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                    cfg_sts.channel = away_lfr_data.label{1};
                    this_p1_sts = ft_spiketriggeredspectrum(cfg_sts, this_p1_data);
                    this_p2_sts = ft_spiketriggeredspectrum(cfg_sts, this_p2_data);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_sts.fourierspctrm{1}),1)) || ...
                            ~isempty(find(isnan(this_p2_sts.fourierspctrm{1}),1))
                        flag_nan_in_split = true;
                        break;
                    end
                    all_freqs(iSplit,:) = this_p1_sts.freq;
                    all_sts_vals(1,iSplit,:) = nanmean(sq(abs(this_p1_sts.fourierspctrm{1})));
                    all_sts_vals(2,iSplit,:) = nanmean(sq(abs(this_p2_sts.fourierspctrm{1})));
                    
                    % Calculate and save PPC
                    cfg_ppc               = [];
                    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg_ppc.spikechannel  = this_sts.label;
                    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                    cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                    cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                    this_p1_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p1_sts);
                    this_p2_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p2_sts);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_ppc.ppc0),1)) || ...
                         ~isempty(find(isnan(this_p2_ppc.ppc0),1))   
                        flag_nan_in_split = true;
                        break;
                    end
                    all_ppc(1,iSplit,:) = this_p1_ppc.ppc0';
                    all_ppc(2,iSplit,:) = this_p2_ppc.ppc0'; 
                end
                if flag_nan_in_split
                   od.msn_res.away_spec{iM}.flag_no_control_split = true;
                   break;
                end
                
                % Save averages and std
                od.msn_res.away_p1_spec{iM}.sta_time = od.msn_res.away_spec{iM}.sta_time;
                od.msn_res.away_p2_spec{iM}.sta_time = od.msn_res.away_spec{iM}.sta_time;
                od.msn_res.away_p1_spec{iM}.freqs = squeeze(mean(all_freqs(1,:,:),2))';
                od.msn_res.away_p2_spec{iM}.freqs = squeeze(mean(all_freqs(2,:,:),2))';
                od.msn_res.away_p1_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(1,:,:),2))';
                od.msn_res.away_p2_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(2,:,:),2))';
                od.msn_res.away_p1_spec{iM}.sd_sta = squeeze(std(all_sta_vals(1,:,:)))';
                od.msn_res.away_p2_spec{iM}.sd_sta = squeeze(std(all_sta_vals(2,:,:)))';
                od.msn_res.away_p1_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(1,:,:),2))';
                od.msn_res.away_p2_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(2,:,:),2))';
                od.msn_res.away_p1_spec{iM}.sd_sts = squeeze(std(all_sts_vals(1,:,:)))';
                od.msn_res.away_p2_spec{iM}.sd_sts = squeeze(std(all_sts_vals(2,:,:)))';
                od.msn_res.away_p1_spec{iM}.mean_ppc = squeeze(mean(all_ppc(1,:,:),2))';
                od.msn_res.away_p2_spec{iM}.mean_ppc = squeeze(mean(all_ppc(2,:,:),2))';
                od.msn_res.away_p1_spec{iM}.sd_ppc = squeeze(std(all_ppc(1,:,:)))';
                od.msn_res.away_p2_spec{iM}.sd_ppc = squeeze(std(all_ppc(2,:,:)))';
                od.msn_res.away_p1_spec{iM}.spk_count = round(mean(all_spk_count(1,:)));
                od.msn_res.away_p2_spec{iM}.spk_count = round(mean(all_spk_count(2,:)));
            else
                od.msn_res.away_spec{iM}.flag_no_control_split = true;
            end
        end
    end
    
    % Get the msn distributions
    od.msn_onTrack_dist = [];
    od.msn_near_dist = [];
    od.msn_near_lfr_dist = [];
    od.msn_near_hfr_dist = [];
    od.msn_near_p1_dist = [];
    od.msn_near_p2_dist = [];
    od.msn_away_dist = [];
    od.msn_away_lfr_dist = [];
    od.msn_away_hfr_dist = [];
    od.msn_away_p1_dist = [];
    od.msn_away_p2_dist = [];
    
    if isfield(od,'msn_res')
        for iM = 1:length(od.msn_res.onTrack_spec)
            od.msn_onTrack_dist = [od.msn_onTrack_dist od.msn_res.onTrack_spec{iM}.spk_count];
            if ~od.msn_res.near_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_spec{iM}.flag_nansts & ~od.msn_res.near_spec{iM}.flag_nanppc
                od.msn_near_dist = [od.msn_near_dist od.msn_res.near_spec{iM}.spk_count];
                if ~od.msn_res.near_lfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_lfr_spec{iM}.flag_nansts & ~od.msn_res.near_lfr_spec{iM}.flag_nanppc
                    od.msn_near_lfr_dist = [od.msn_near_lfr_dist od.msn_res.near_lfr_spec{iM}.spk_count];
                end
                if ~od.msn_res.near_hfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.near_hfr_spec{iM}.flag_nansts & ~od.msn_res.near_hfr_spec{iM}.flag_nanppc
                    od.msn_near_hfr_dist = [od.msn_near_hfr_dist od.msn_res.near_hfr_spec{iM}.spk_count];
                end
                if ~od.msn_res.near_spec{iM}.flag_no_control_split
                    od.msn_near_p1_dist = [od.msn_near_p1_dist od.msn_res.near_p1_spec{iM}.spk_count];
                    od.msn_near_p2_dist = [od.msn_near_p2_dist od.msn_res.near_p2_spec{iM}.spk_count];
                end
            end
            if ~od.msn_res.away_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_spec{iM}.flag_nansts & ~od.msn_res.away_spec{iM}.flag_nanppc
                od.msn_away_dist = [od.msn_away_dist od.msn_res.away_spec{iM}.spk_count];
                if ~od.msn_res.away_lfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_lfr_spec{iM}.flag_nansts & ~od.msn_res.away_lfr_spec{iM}.flag_nanppc
                    od.msn_away_lfr_dist = [od.msn_away_lfr_dist od.msn_res.away_lfr_spec{iM}.spk_count];
                end
                if ~od.msn_res.away_hfr_spec{iM}.flag_tooFewSpikes & ~od.msn_res.away_hfr_spec{iM}.flag_nansts & ~od.msn_res.away_hfr_spec{iM}.flag_nanppc
                    od.msn_away_hfr_dist = [od.msn_away_hfr_dist od.msn_res.away_hfr_spec{iM}.spk_count];
                end
               if ~od.msn_res.away_spec{iM}.flag_no_control_split
                    od.msn_away_p1_dist = [od.msn_away_p1_dist od.msn_res.away_p1_spec{iM}.spk_count];
                    od.msn_away_p2_dist = [od.msn_away_p2_dist od.msn_res.away_p2_spec{iM}.spk_count];
                end
            end
        end
    end
    
    % Calculate spectral measures for all FSIs
    all_fsi = find(od.cell_type == 2);
    for iM = 1:length(all_fsi)
        iC  = all_fsi(iM);
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
        od.fsi_res.onTrack_spec{iM}.spk_count = sum(on_track_data.trial{1}(2,:));
        this_sta = ft_spiketriggeredaverage(cfg_ft, on_track_data);
        od.fsi_res.onTrack_spec{iM}.sta_time = this_sta.time;
        od.fsi_res.onTrack_spec{iM}.sta_vals = this_sta.avg(:,:)';
        od.fsi_res.onTrack_spec{iM}.flag_unequalSpikes = this_flag;
        
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
        od.fsi_res.onTrack_spec{iM}.freqs = this_sts.freq;
        od.fsi_res.onTrack_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
        od.fsi_res.onTrack_spec{iM}.flag_nansts = this_flag;
        
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
        od.fsi_res.onTrack_spec{iM}.ppc = this_ppc.ppc0';
        od.fsi_res.onTrack_spec{iM}.flag_nanppc = this_flag;

        % Calculate and save subsampled measures for onTrack
        % When an FSI has fewer spikes than at least one MSN, do NOT
        % subsample and include all the spikes
        if isempty(od.msn_onTrack_dist)
            od.fsi_res.onTrack_spec{iM}.flag_no_subsampling = true;
        elseif od.fsi_res.onTrack_spec{iM}.spk_count < max(od.msn_onTrack_dist)
            od.fsi_res.onTrack_spec{iM}.flag_no_subsampling = true;
        elseif isfield(od.fsi_res.onTrack_spec{iM},'flag_nansts') && od.fsi_res.onTrack_spec{iM}.flag_nansts
            od.fsi_res.onTrack_spec{iM}.flag_no_subsampling = true;
        elseif isfield(od.fsi_res.onTrack_spec{iM},'flag_nanppc') && od.fsi_res.onTrack_spec{iM}.flag_nanppc
            od.fsi_res.onTrack_spec{iM}.flag_no_subsampling = true;
        else
            od.fsi_res.onTrack_spec{iM}.flag_no_subsampling = false;
            temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.onTrack_spec{iM}.sts_vals));
            temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.onTrack_spec{iM}.ppc));
            for iS = 1:cfg_master.num_subsamples
                s_factor = 1/length(od.msn_onTrack_dist);
                choice = floor(rand()/s_factor)+1;
                sub_idx = randsample(1:od.fsi_res.onTrack_spec{iM}.spk_count, od.msn_onTrack_dist(choice));
                sub_idx = sort(sub_idx); sub_sts = this_sts;
                sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                temp_ppc_vals(iS,:) = sub_ppc.ppc0';
            end
            od.fsi_res.onTrack_spec{iM}.subsampled_sts = mean(temp_sts_vals,1);
            od.fsi_res.onTrack_spec{iM}.subsampled_ppc = mean(temp_ppc_vals,1);
        end
        
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
            od.fsi_res.near_spec{iM}.flag_tooFewSpikes = true;
        else
            od.fsi_res.near_spec{iM}.spk_count = spk_count2;
            od.fsi_res.near_spec{iM}.flag_tooFewSpikes = false;
            if spk_count1 ~= spk_count2
                this_flag = true;
                warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in near-Reward trials', ...
                    spk_count1, spk_count2)

            end
            % Calculate and save STA
            this_sta = ft_spiketriggeredaverage(cfg_ft, near_data);
            od.fsi_res.near_spec{iM}.sta_time = this_sta.time;
            od.fsi_res.near_spec{iM}.sta_vals = this_sta.avg(:,:)';
            od.fsi_res.near_spec{iM}.flag_unequalSpikes = this_flag;

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
            od.fsi_res.near_spec{iM}.freqs = this_sts.freq;
            od.fsi_res.near_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
            od.fsi_res.near_spec{iM}.flag_nansts = this_flag;

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
            od.fsi_res.near_spec{iM}.ppc = this_ppc.ppc0';
            od.fsi_res.near_spec{iM}.flag_nanppc = this_flag;
            
            % Calculate and save subsampled measures for near Reward trials
            % When an FSI has fewer spikes than at least one MSN, do NOT
            % subsample and include all the spikes
            if isempty(od.msn_near_dist)
                od.fsi_res.near_spec{iM}.flag_no_subsampling = true;
            elseif od.fsi_res.near_spec{iM}.spk_count < max(od.msn_near_dist)
                od.fsi_res.near_spec{iM}.flag_no_subsampling = true;
            elseif isfield(od.fsi_res.near_spec{iM},'flag_nansts') && od.fsi_res.near_spec{iM}.flag_nansts
                od.fsi_res.near_spec{iM}.flag_no_subsampling = true;
            elseif isfield(od.fsi_res.near_spec{iM},'flag_nanppc') && od.fsi_res.near_spec{iM}.flag_nanppc
                od.fsi_res.near_spec{iM}.flag_no_subsampling = true;
            elseif isfield(od.fsi_res.onTrack_spec{iM},'flag_no_subsampling') && od.fsi_res.onTrack_spec{iM}.flag_no_subsampling
                od.fsi_res.near_spec{iM}.flag_no_subsampling = true;
            else
                od.fsi_res.near_spec{iM}.flag_no_subsampling = false;
                temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_spec{iM}.sts_vals));
                temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_spec{iM}.ppc));
                for iS = 1:cfg_master.num_subsamples
                    s_factor = 1/length(od.msn_near_dist);
                    choice = floor(rand()/s_factor)+1;
                    sub_idx = randsample(1:od.fsi_res.near_spec{iM}.spk_count, od.msn_near_dist(choice));
                    sub_idx = sort(sub_idx); sub_sts = this_sts;
                    sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                    sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                    sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                    temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                    sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                    temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                end
                od.fsi_res.near_spec{iM}.subsampled_sts = mean(temp_sts_vals,1);
                od.fsi_res.near_spec{iM}.subsampled_ppc = mean(temp_ppc_vals,1);
            end
            
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
            % Find out the firing rate threshold to split the trials such that
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
            
            % Divide trials into nControlSplits partitions with nearly equal spikes but not on
            % the basis of fr_threshold. 
            % Algo used: Use randomly generated splits and accept a split only if
            % 1) difference between splits is <= 2x the diffference of the
            % FR based split
            % 2) not the same as an experimental split
            % 3) not the same as an already found split
            % Wait for 100,000,000 iterations, if you haven't found
            % nControlSplits valid splits by then, take what you have!
            
            valid_splits  = false(cfg_master.nControlSplits,tcount);
            all_tspikes_mat = cell2mat(all_tspikes);
            last_valid_split = 0;
            for iRand = 1:100000000 % wait till 100 million iterations
                A = 1:tcount;
                ndiv = 2;
                this_idx = sort([1 randperm(length(A)-1, ndiv-1)+1 length(A)+1]);
                for k1 = 1:length(this_idx)-1
                    R{k1} = A(this_idx(k1):this_idx(k1+1)-1);
                end  
                this_perm = randperm(tcount);
                this_split = false(1,tcount);
                this_split(this_perm(R{1})) = true;
                this_split_dif = abs(sum(all_tspikes_mat(this_split)) - sum(all_tspikes_mat(~this_split)));
                % reject split of split dif is grater than threshold
                if this_split_dif > dif_min
                   continue; 
                end
                % reject if split exactly the same as hypothesis split
                if sum(this_split == hfr_trials) == tcount | sum(this_split == lfr_trials) == tcount
                   continue;
                end
                % reject if split the same as previously found valid split
                flag_repeat_split = false;
                for iCheck = 1:1:last_valid_split
                    if sum(this_split == valid_splits(iCheck)) == tcount | sum(~this_split == valid_splits(iCheck)) == tcount
                        flag_repeat_spilt = true;
                        break;
                    end
                end
                if flag_repeat_split
                    continue;
                end
                % If you have made it till here, you found a valid_split!
                last_valid_split = last_valid_split + 1;
                valid_splits(last_valid_split,:) = this_split;
                % If 100 valid splits are found, get out of this!
                if last_valid_split == cfg_master.nControlSplits
                    break;
                end
            end
              
            od.fsi_res.near_spec{iM}.mfr = mfr;
            od.fsi_res.near_spec{iM}.fr_thresh = fr_thresh;
            od.fsi_res.near_spec{iM}.trial_spk_count = cell2mat(all_tspikes);
            od.fsi_res.near_spec{iM}.valid_split_count = last_valid_split;
            od.fsi_res.near_spec{iM}.valid_splits = valid_splits;
            od.fsi_res.near_spec{iM}.randIters = iRand;

            % Calculate and Save all spec results for near_hfr_data
            cfg_near_hfr_trials.trl = cfg_near_trials.trl(hfr_trials,:);
            near_hfr_data = ft_redefinetrial(cfg_near_hfr_trials, this_data);
            
            % Skip if no spks in near_hfr
            spk_count = 0;
            for iT = 1:length(near_hfr_data.trial)
                spk_count = spk_count + sum(near_hfr_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes2
                od.fsi_res.near_hfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.fsi_res.near_hfr_spec{iM}.spk_count = spk_count;
                od.fsi_res.near_hfr_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, near_hfr_data);
                od.fsi_res.near_hfr_spec{iM}.sta_time = this_sta.time;
                od.fsi_res.near_hfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

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
                od.fsi_res.near_hfr_spec{iM}.flag_nansts = this_flag;
                od.fsi_res.near_hfr_spec{iM}.freqs = this_sts.freq;
                od.fsi_res.near_hfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));

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
                od.fsi_res.near_hfr_spec{iM}.ppc = this_ppc.ppc0';
                od.fsi_res.near_hfr_spec{iM}.flag_nanppc = this_flag;
                
                % Calculate and save subsampled measures for near Reward hfr trials
                % When an FSI has fewer spikes than at least one MSN, do NOT
                % subsample and include all the spikes
                if isempty(od.msn_near_hfr_dist)
                    od.fsi_res.near_hfr_spec{iM}.flag_no_subsampling = true;
                elseif od.fsi_res.near_hfr_spec{iM}.spk_count < max(od.msn_near_hfr_dist)
                    od.fsi_res.near_hfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.near_hfr_spec{iM},'flag_nansts') && od.fsi_res.near_hfr_spec{iM}.flag_nansts
                    od.fsi_res.near_hfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.near_hfr_spec{iM},'flag_nanppc') && od.fsi_res.near_hfr_spec{iM}.flag_nanppc
                    od.fsi_res.near_hfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.near_spec{iM},'flag_no_subsampling') && od.fsi_res.near_spec{iM}.flag_no_subsampling
                    od.fsi_res.near_hfr_spec{iM}.flag_no_subsampling = true;
                else
                    od.fsi_res.near_hfr_spec{iM}.flag_no_subsampling = false;
                    temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_hfr_spec{iM}.sts_vals));
                    temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_hfr_spec{iM}.ppc));
                    for iS = 1:cfg_master.num_subsamples
                        s_factor = 1/length(od.msn_near_hfr_dist);
                        choice = floor(rand()/s_factor)+1;
                        sub_idx = randsample(1:od.fsi_res.near_hfr_spec{iM}.spk_count, od.msn_near_hfr_dist(choice));
                        sub_idx = sort(sub_idx); sub_sts = this_sts;
                        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                        temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                        temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                    end
                    od.fsi_res.near_hfr_spec{iM}.subsampled_sts = mean(temp_sts_vals,1);
                    od.fsi_res.near_hfr_spec{iM}.subsampled_ppc = mean(temp_ppc_vals,1);
                end
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
                od.fsi_res.near_lfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.fsi_res.near_lfr_spec{iM}.spk_count = spk_count;
                od.fsi_res.near_lfr_spec{iM}.flag_tooFewSpikes = false;
                
                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, near_lfr_data);
                od.fsi_res.near_lfr_spec{iM}.sta_time = this_sta.time;
                od.fsi_res.near_lfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

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
                od.fsi_res.near_lfr_spec{iM}.freqs = this_sts.freq;
                od.fsi_res.near_lfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                od.fsi_res.near_lfr_spec{iM}.flag_nansts = this_flag;

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
                od.fsi_res.near_lfr_spec{iM}.ppc = this_ppc.ppc0';
                od.fsi_res.near_lfr_spec{iM}.flag_nanppc = this_flag;
                
                % Calculate and save subsampled measures for near Reward lfr trials
                % When an FSI has fewer spikes than at least one MSN, do NOT
                % subsample and include all the spikes
                if isempty(od.msn_near_lfr_dist)
                    od.fsi_res.near_lfr_spec{iM}.flag_no_subsampling = true;
                elseif od.fsi_res.near_lfr_spec{iM}.spk_count < max(od.msn_near_lfr_dist)
                    od.fsi_res.near_lfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.near_lfr_spec{iM},'flag_nansts') && od.fsi_res.near_lfr_spec{iM}.flag_nansts
                    od.fsi_res.near_lfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.near_lfr_spec{iM},'flag_nanppc') && od.fsi_res.near_lfr_spec{iM}.flag_nanppc
                    od.fsi_res.near_lfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.near_spec{iM},'flag_no_subsampling') && od.fsi_res.near_spec{iM}.flag_no_subsampling
                    od.fsi_res.near_lfr_spec{iM}.flag_no_subsampling = true;
                else
                    od.fsi_res.near_lfr_spec{iM}.flag_no_subsampling = false;
                    temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_lfr_spec{iM}.sts_vals));
                    temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_lfr_spec{iM}.ppc));
                    for iS = 1:cfg_master.num_subsamples
                        s_factor = 1/length(od.msn_near_lfr_dist);
                        choice = floor(rand()/s_factor)+1;
                        sub_idx = randsample(1:od.fsi_res.near_lfr_spec{iM}.spk_count, od.msn_near_lfr_dist(choice));
                        sub_idx = sort(sub_idx); sub_sts = this_sts;
                        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                        temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                        temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                    end
                    od.fsi_res.near_lfr_spec{iM}.subsampled_sts = mean(temp_sts_vals,1);
                    od.fsi_res.near_lfr_spec{iM}.subsampled_ppc = mean(temp_ppc_vals,1);
                end
            end
            
            % Do control split stuff only if lfr/mfr splits are perfect
            % i.e., enough spikes and no nan flags
            
            if ~od.fsi_res.near_lfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.fsi_res.near_hfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.fsi_res.near_lfr_spec{iM}.flag_nansts && ...
                    ~od.fsi_res.near_lfr_spec{iM}.flag_nanppc  && ...
                    ~od.fsi_res.near_hfr_spec{iM}.flag_nansts && ...
                    ~od.fsi_res.near_hfr_spec{iM}.flag_nanppc && ...
                    od.fsi_res.near_spec{iM}.valid_split_count == cfg_master.nControlSplits
                
                od.fsi_res.near_spec{iM}.flag_no_control_split = false;
                
                % Create temporary arrays to average later
                all_sta_vals = zeros([2, cfg_master.nControlSplits, size(od.fsi_res.near_spec{iM}.sta_vals)]);
                all_freqs = zeros([cfg_master.nControlSplits, size(od.fsi_res.near_spec{iM}.freqs)]);
                all_sts_vals = zeros([2, cfg_master.nControlSplits, length(od.fsi_res.near_spec{iM}.sts_vals)]);
                all_ppc = zeros([2, cfg_master.nControlSplits, size(od.fsi_res.near_spec{iM}.ppc)]);
                all_spk_count = zeros(2, cfg_master.nControlSplits);
                all_subsampled_sts = zeros([2, cfg_master.nControlSplits, length(od.fsi_res.near_spec{iM}.sts_vals)]);
                all_subsampled_ppc = zeros([2, cfg_master.nControlSplits, size(od.fsi_res.near_spec{iM}.ppc)]);

                % For each split, fill in the tables
                % If any one of them is nan, exit loop, set
                % flag_no_control_split as true and move on to the next
                % cell    
                flag_nan_in_split = false;      
                for iSplit = 1:cfg_master.nControlSplits
  
                    % Calculate and Save all spec results for near_lfr_data
                    this_p1_trials = find(od.fsi_res.near_spec{iM}.valid_splits(iSplit,:));
                    this_p1_cfg.trl = cfg_near_trials.trl(this_p1_trials,:);
                    this_p1_data = ft_redefinetrial(this_p1_cfg, this_data);
                    this_p2_trials = find(~od.fsi_res.near_spec{iM}.valid_splits(iSplit,:));
                    this_p2_cfg.trl = cfg_near_trials.trl(this_p2_trials,:);
                    this_p2_data = ft_redefinetrial(this_p2_cfg, this_data);
                    
                    p1_spk_count = 0;
                    for iT = 1:length(this_p1_data.trial)
                        p1_spk_count = p1_spk_count + sum(this_p1_data.trial{iT}(2,:)); 
                    end
                    
                    p2_spk_count = 0;
                    for iT = 1:length(this_p2_data.trial)
                        p2_spk_count = p2_spk_count + sum(this_p2_data.trial{iT}(2,:)); 
                    end
                    
                    all_spk_count(1,iSplit) = p1_spk_count;
                    all_spk_count(2,iSplit) = p2_spk_count;

                    % Calculate and save STA
                    this_p1_sta = ft_spiketriggeredaverage(cfg_ft, this_p1_data);
                    this_p2_sta = ft_spiketriggeredaverage(cfg_ft, this_p2_data);
                    all_sta_vals(1, iSplit, :) = this_p1_sta.avg(:,:)';
                    all_sta_vals(2, iSplit, :) = this_p2_sta.avg(:,:)';
                    
                    % Calculate and save STS
                    cfg_sts.method = 'mtmconvol';
                    cfg_sts.foi = 1:1:100;
                    cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                    cfg_sts.taper = 'hanning';
                    cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                    cfg_sts.channel = near_lfr_data.label{1};
                    this_p1_sts = ft_spiketriggeredspectrum(cfg_sts, this_p1_data);
                    this_p2_sts = ft_spiketriggeredspectrum(cfg_sts, this_p2_data);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_sts.fourierspctrm{1}),1)) || ...
                            ~isempty(find(isnan(this_p2_sts.fourierspctrm{1}),1))
                        flag_nan_in_split = true;
                        break;
                    end
                    all_freqs(iSplit,:) = this_p1_sts.freq;
                    all_sts_vals(1,iSplit,:) = nanmean(sq(abs(this_p1_sts.fourierspctrm{1})));
                    all_sts_vals(2,iSplit,:) = nanmean(sq(abs(this_p2_sts.fourierspctrm{1})));
                    
                    % Calculate and save PPC
                    cfg_ppc               = [];
                    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg_ppc.spikechannel  = this_sts.label;
                    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                    cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                    cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                    this_p1_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p1_sts);
                    this_p2_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p2_sts);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_ppc.ppc0),1)) || ...
                         ~isempty(find(isnan(this_p2_ppc.ppc0),1))   
                        flag_nan_in_split = true;
                        break;
                    end
                    all_ppc(1,iSplit,:) = this_p1_ppc.ppc0';
                    all_ppc(2,iSplit,:) = this_p2_ppc.ppc0'; 

					% Calculate and save subsampled measures for near Reward p1 trials
					% When an FSI has fewer spikes than at least one MSN, do NOT
					% subsample and include all the spikes
					if isempty(od.msn_near_p1_dist)
						od.fsi_res.near_p1_spec{iM}.flag_no_subsampling = true;
                    elseif p1_spk_count < max(od.msn_near_p1_dist)
                        od.fsi_res.near_p1_spec{iM}.flag_no_subsampling = true;
                    elseif isfield(od.fsi_res.near_spec{iM},'flag_no_subsampling') && od.fsi_res.near_spec{iM}.flag_no_subsampling
                        od.fsi_res.near_p1_spec{iM}.flag_no_subsampling = true;
					else
						od.fsi_res.near_p1_spec{iM}.flag_no_subsampling = false;
 							temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_spec{iM}.sts_vals));
							temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_spec{iM}.ppc));
						for iS = 1:cfg_master.num_subsamples
						    s_factor = 1/length(od.msn_near_p1_dist);
						    choice = floor(rand()/s_factor)+1;
						    sub_idx = randsample(1:p1_spk_count, od.msn_near_p1_dist(choice));
						    sub_idx = sort(sub_idx); 
						    sub_sts = this_p1_sts;
						    sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
						    sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
						    sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
						    temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
						    sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
						    temp_ppc_vals(iS,:) = sub_ppc.ppc0';
						end
						all_subsampled_sts(1,iSplit,:) = mean(temp_sts_vals,1);
						all_subsampled_ppc(1,iSplit,:) = mean(temp_ppc_vals,1);
					end

					% Calculate and save subsampled measures for near Reward p2 trials
                    % When an FSI has fewer spikes than at least one MSN, do NOT
                    % subsample and include all the spikes
                    if isempty(od.msn_near_p2_dist)
						od.fsi_res.near_p2_spec{iM}.flag_no_subsampling = true;
                    elseif p2_spk_count < max(od.msn_near_p2_dist)
                        od.fsi_res.near_p2_spec{iM}.flag_no_subsampling = true;
                    elseif isfield(od.fsi_res.near_spec{iM},'flag_no_subsampling') && od.fsi_res.near_spec{iM}.flag_no_subsampling
                        od.fsi_res.near_p2_spec{iM}.flag_no_subsampling = true;
					else
                        od.fsi_res.near_p2_spec{iM}.flag_no_subsampling = false;
                            temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_spec{iM}.sts_vals));
                            temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.near_spec{iM}.ppc));
                        for iS = 1:cfg_master.num_subsamples
                            s_factor = 1/length(od.msn_near_p2_dist);
                            choice = floor(rand()/s_factor)+1;
                            sub_idx = randsample(1:p2_spk_count, od.msn_near_p2_dist(choice));
                            sub_idx = sort(sub_idx); 
                            sub_sts = this_p2_sts;
                            sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                            sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                            sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                            temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                            sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                            temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                        end
                        all_subsampled_sts(2,iSplit,:) = mean(temp_sts_vals,1);
                        all_subsampled_ppc(2,iSplit,:) = mean(temp_ppc_vals,1);
                    end
                end
                if flag_nan_in_split
                   od.fsi_res.near_spec{iM}.flag_no_control_split = true;
                   break;
                end
                % Save averages and std
                od.fsi_res.near_p1_spec{iM}.sta_time = od.fsi_res.near_spec{iM}.sta_time;
                od.fsi_res.near_p2_spec{iM}.sta_time = od.fsi_res.near_spec{iM}.sta_time;
                od.fsi_res.near_p1_spec{iM}.freqs = squeeze(mean(all_freqs(1,:,:),2))';
                od.fsi_res.near_p2_spec{iM}.freqs = squeeze(mean(all_freqs(2,:,:),2))';
                od.fsi_res.near_p1_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(1,:,:),2))';
                od.fsi_res.near_p2_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(2,:,:),2))';
                od.fsi_res.near_p1_spec{iM}.sd_sta = squeeze(std(all_sta_vals(1,:,:)))';
                od.fsi_res.near_p2_spec{iM}.sd_sta = squeeze(std(all_sta_vals(2,:,:)))';
                od.fsi_res.near_p1_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(1,:,:),2))';
                od.fsi_res.near_p2_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(2,:,:),2))';
                od.fsi_res.near_p1_spec{iM}.sd_sts = squeeze(std(all_sts_vals(1,:,:)))';
                od.fsi_res.near_p2_spec{iM}.sd_sts = squeeze(std(all_sts_vals(2,:,:)))';
                od.fsi_res.near_p1_spec{iM}.mean_ppc = squeeze(mean(all_ppc(1,:,:),2))';
                od.fsi_res.near_p2_spec{iM}.mean_ppc = squeeze(mean(all_ppc(2,:,:),2))';
                od.fsi_res.near_p1_spec{iM}.sd_ppc = squeeze(std(all_ppc(1,:,:)))';
                od.fsi_res.near_p2_spec{iM}.sd_ppc = squeeze(std(all_ppc(2,:,:)))';
                od.fsi_res.near_p1_spec{iM}.mean_subsampled_sts = squeeze(mean(all_subsampled_sts(1,:,:),2))';
                od.fsi_res.near_p2_spec{iM}.mean_subsampled_sts = squeeze(mean(all_subsampled_sts(2,:,:),2))';
                od.fsi_res.near_p1_spec{iM}.sd_subsampled_sts = squeeze(std(all_subsampled_sts(1,:,:)))';
                od.fsi_res.near_p2_spec{iM}.sd_subsampled_sts = squeeze(std(all_subsampled_sts(2,:,:)))';
                od.fsi_res.near_p1_spec{iM}.mean_subsampled_ppc = squeeze(mean(all_subsampled_ppc(1,:,:),2))';
                od.fsi_res.near_p2_spec{iM}.mean_subsampled_ppc = squeeze(mean(all_subsampled_ppc(2,:,:),2))';
                od.fsi_res.near_p1_spec{iM}.sd_subsampled_ppc = squeeze(std(all_subsampled_ppc(1,:,:)))';
                od.fsi_res.near_p2_spec{iM}.sd_subsampled_ppc = squeeze(std(all_subsampled_ppc(2,:,:)))';
                od.fsi_res.near_p1_spec{iM}.spk_count = round(mean(all_spk_count(1,:)));
                od.fsi_res.near_p2_spec{iM}.spk_count = round(mean(all_spk_count(2,:)));
            else
                od.fsi_res.near_spec{iM}.flag_no_control_split = true;
            end
        end 
        % For away reward_trials  
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
        if spk_count2 <  cfg_master.nMinSpikes2
            od.fsi_res.away_spec{iM}.flag_tooFewSpikes = true;
        else
            od.fsi_res.away_spec{iM}.spk_count = spk_count2;
            od.fsi_res.away_spec{iM}.flag_tooFewSpikes = false;
            if spk_count1 ~= spk_count2
                this_flag = true;
                warning('ft_redefinetrial has %d spikes but restrict shows %d spikes in away-Reward trials', ...
                    spk_count1, spk_count2)

            end
            % Calculate and save STA
            this_sta = ft_spiketriggeredaverage(cfg_ft, away_data);
            od.fsi_res.away_spec{iM}.sta_time = this_sta.time;
            od.fsi_res.away_spec{iM}.sta_vals = this_sta.avg(:,:)';
            od.fsi_res.away_spec{iM}.flag_unequalSpikes = this_flag;

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
            od.fsi_res.away_spec{iM}.freqs = this_sts.freq;
            od.fsi_res.away_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
            od.fsi_res.away_spec{iM}.flag_nansts = this_flag;

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
            od.fsi_res.away_spec{iM}.ppc = this_ppc.ppc0';
            od.fsi_res.away_spec{iM}.flag_nanppc = this_flag;

            % Calculate and save subsampled measures for away Reward trials
            % When an FSI has fewer spikes than at least one MSN, do NOT
            % subsample and include all the spikes
            if isempty(od.msn_away_dist)
                od.fsi_res.away_spec{iM}.flag_no_subsampling = true;
            elseif od.fsi_res.away_spec{iM}.spk_count < max(od.msn_away_dist)
                od.fsi_res.away_spec{iM}.flag_no_subsampling = true;
            elseif isfield(od.fsi_res.away_spec{iM},'flag_nansts') && od.fsi_res.away_spec{iM}.flag_nansts
                od.fsi_res.away_spec{iM}.flag_no_subsampling = true;
            elseif isfield(od.fsi_res.away_spec{iM},'flag_nanppc') && od.fsi_res.away_spec{iM}.flag_nanppc
                od.fsi_res.away_spec{iM}.flag_no_subsampling = true;
            elseif isfield(od.fsi_res.onTrack_spec{iM},'flag_no_subsampling') && od.fsi_res.onTrack_spec{iM}.flag_no_subsampling
                od.fsi_res.away_spec{iM}.flag_no_subsampling = true;
            else
                od.fsi_res.away_spec{iM}.flag_no_subsampling = false;
                temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_spec{iM}.sts_vals));
                temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_spec{iM}.ppc));
                for iS = 1:cfg_master.num_subsamples
                    s_factor = 1/length(od.msn_away_dist);
                    choice = floor(rand()/s_factor)+1;
                    sub_idx = randsample(1:od.fsi_res.away_spec{iM}.spk_count, od.msn_away_dist(choice));
                    sub_idx = sort(sub_idx); sub_sts = this_sts;
                    sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                    sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                    sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                    temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                    sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                    temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                end
                od.fsi_res.away_spec{iM}.subsampled_sts = mean(temp_sts_vals,1);
                od.fsi_res.away_spec{iM}.subsampled_ppc = mean(temp_ppc_vals,1);
            end

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
            % Find out the firing rate threshold to split the trials such that
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


           % Divide trials into nControlSplits partitions with awayly equal spikes but not on
            % the basis of fr_threshold. 
            % Algo used: Use randomly generated splits and accept a split only if
            % 1) difference between splits is <= 2x the diffference of the
            % FR based split
            % 2) not the same as an experimental split
            % 3) not the same as an already found split
            % Wait for 100,000,000 iterations, if you haven't found
            % nControlSplits valid splits by then, take what you have!

            valid_splits  = false(cfg_master.nControlSplits,tcount);
            all_tspikes_mat = cell2mat(all_tspikes);
            last_valid_split = 0;
            for iRand = 1:100000000 % wait till 100 million iterations
                A = 1:tcount;
                ndiv = 2;
                this_idx = sort([1 randperm(length(A)-1, ndiv-1)+1 length(A)+1]);
                for k1 = 1:length(this_idx)-1
                    R{k1} = A(this_idx(k1):this_idx(k1+1)-1);
                end  
                this_perm = randperm(tcount);
                this_split = false(1,tcount);
                this_split(this_perm(R{1})) = true;
                this_split_dif = abs(sum(all_tspikes_mat(this_split)) - sum(all_tspikes_mat(~this_split)));
                % reject split of split dif is grater than threshold
                if this_split_dif > dif_min
                   continue; 
                end
                % reject if split exactly the same as hypothesis split
                if sum(this_split == hfr_trials) == tcount | sum(this_split == lfr_trials) == tcount
                   continue;
                end
                % reject if split the same as previously found valid split
                flag_repeat_split = false;
                for iCheck = 1:1:last_valid_split
                    if sum(this_split == valid_splits(iCheck)) == tcount | sum(~this_split == valid_splits(iCheck)) == tcount
                        flag_repeat_spilt = true;
                        break;
                    end
                end
                if flag_repeat_split
                    continue;
                end
                % If you have made it till here, you found a valid_split!
                last_valid_split = last_valid_split + 1;
                valid_splits(last_valid_split,:) = this_split;
                % If 100 valid splits are found, get out of this!
                if last_valid_split == cfg_master.nControlSplits
                    break;
                end
            end

            od.fsi_res.away_spec{iM}.mfr = mfr;
            od.fsi_res.away_spec{iM}.fr_thresh = fr_thresh;
            od.fsi_res.away_spec{iM}.trial_spk_count = cell2mat(all_tspikes);
            od.fsi_res.away_spec{iM}.valid_split_count = last_valid_split;
            od.fsi_res.away_spec{iM}.valid_splits = valid_splits;
            od.fsi_res.away_spec{iM}.randIters = iRand;

            % Calculate and Save all spec results for away_hfr_data
            cfg_away_hfr_trials.trl = cfg_away_trials.trl(hfr_trials,:);
            away_hfr_data = ft_redefinetrial(cfg_away_hfr_trials, this_data);

            % Skip if no spks in away_hfr
            spk_count = 0;
            for iT = 1:length(away_hfr_data.trial)
                spk_count = spk_count + sum(away_hfr_data.trial{iT}(2,:)); 
            end
            if spk_count < cfg_master.nMinSpikes2
                od.fsi_res.away_hfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.fsi_res.away_hfr_spec{iM}.spk_count = spk_count;
                od.fsi_res.away_hfr_spec{iM}.flag_tooFewSpikes = false;

                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, away_hfr_data);
                od.fsi_res.away_hfr_spec{iM}.sta_time = this_sta.time;
                od.fsi_res.away_hfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

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
                od.fsi_res.away_hfr_spec{iM}.flag_nansts = this_flag;
                od.fsi_res.away_hfr_spec{iM}.freqs = this_sts.freq;
                od.fsi_res.away_hfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));

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
                od.fsi_res.away_hfr_spec{iM}.ppc = this_ppc.ppc0';
                od.fsi_res.away_hfr_spec{iM}.flag_nanppc = this_flag;

                % Calculate and save subsampled measures for away Reward hfr trials
                % When an FSI has fewer spikes than at least one MSN, do NOT
                % subsample and include all the spikes
                if isempty(od.msn_away_hfr_dist)
                    od.fsi_res.away_hfr_spec{iM}.flag_no_subsampling = true;
                elseif od.fsi_res.away_hfr_spec{iM}.spk_count < max(od.msn_away_hfr_dist)
                    od.fsi_res.away_hfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.away_hfr_spec{iM},'flag_nansts') && od.fsi_res.away_hfr_spec{iM}.flag_nansts
                    od.fsi_res.away_hfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.away_hfr_spec{iM},'flag_nanppc') && od.fsi_res.away_hfr_spec{iM}.flag_nanppc
                    od.fsi_res.away_hfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.away_spec{iM},'flag_no_subsampling') && od.fsi_res.away_spec{iM}.flag_no_subsampling
                    od.fsi_res.away_hfr_spec{iM}.flag_no_subsampling = true;
                else
                    od.fsi_res.away_hfr_spec{iM}.flag_no_subsampling = false;
                    temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_hfr_spec{iM}.sts_vals));
                    temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_hfr_spec{iM}.ppc));
                    for iS = 1:cfg_master.num_subsamples
                        s_factor = 1/length(od.msn_away_hfr_dist);
                        choice = floor(rand()/s_factor)+1;
                        sub_idx = randsample(1:od.fsi_res.away_hfr_spec{iM}.spk_count, od.msn_away_hfr_dist(choice));
                        sub_idx = sort(sub_idx); sub_sts = this_sts;
                        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                        temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                        temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                    end
                    od.fsi_res.away_hfr_spec{iM}.subsampled_sts = mean(temp_sts_vals,1);
                    od.fsi_res.away_hfr_spec{iM}.subsampled_ppc = mean(temp_ppc_vals,1);
                end
            end

            % Calculate and Save all spec results for away_lfr_data
            cfg_away_lfr_trials.trl = cfg_away_trials.trl(lfr_trials,:);
            away_lfr_data = ft_redefinetrial(cfg_away_lfr_trials, this_data);

            % Skip if no spks in away_lfr
            spk_count = 0;
            for iT = 1:length(away_lfr_data.trial)
                spk_count = spk_count + sum(away_lfr_data.trial{iT}(2,:)); 
            end
            if spk_count <  cfg_master.nMinSpikes3
                od.fsi_res.away_lfr_spec{iM}.flag_tooFewSpikes = true;
            else
                od.fsi_res.away_lfr_spec{iM}.spk_count = spk_count;
                od.fsi_res.away_lfr_spec{iM}.flag_tooFewSpikes = false;

                % Calculate and save STA
                this_sta = ft_spiketriggeredaverage(cfg_ft, away_lfr_data);
                od.fsi_res.away_lfr_spec{iM}.sta_time = this_sta.time;
                od.fsi_res.away_lfr_spec{iM}.sta_vals = this_sta.avg(:,:)';

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
                od.fsi_res.away_lfr_spec{iM}.freqs = this_sts.freq;
                od.fsi_res.away_lfr_spec{iM}.sts_vals = nanmean(sq(abs(this_sts.fourierspctrm{1})));
                od.fsi_res.away_lfr_spec{iM}.flag_nansts = this_flag;

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
                od.fsi_res.away_lfr_spec{iM}.ppc = this_ppc.ppc0';
                od.fsi_res.away_lfr_spec{iM}.flag_nanppc = this_flag;

                % Calculate and save subsampled measures for away Reward lfr trials
                % When an FSI has fewer spikes than at least one MSN, do NOT
                % subsample and include all the spikes
                if isempty(od.msn_away_lfr_dist)
                    od.fsi_res.away_lfr_spec{iM}.flag_no_subsampling = true;
                elseif od.fsi_res.away_lfr_spec{iM}.spk_count < max(od.msn_away_lfr_dist)
                    od.fsi_res.away_lfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.away_lfr_spec{iM},'flag_nansts') && od.fsi_res.away_lfr_spec{iM}.flag_nansts
                    od.fsi_res.away_lfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.away_lfr_spec{iM},'flag_nanppc') && od.fsi_res.away_lfr_spec{iM}.flag_nanppc
                    od.fsi_res.away_lfr_spec{iM}.flag_no_subsampling = true;
                elseif isfield(od.fsi_res.away_spec{iM},'flag_no_subsampling') && od.fsi_res.away_spec{iM}.flag_no_subsampling
                    od.fsi_res.away_lfr_spec{iM}.flag_no_subsampling = true;
                else
                    od.fsi_res.away_lfr_spec{iM}.flag_no_subsampling = false;
                    temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_lfr_spec{iM}.sts_vals));
                    temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_lfr_spec{iM}.ppc));
                    for iS = 1:cfg_master.num_subsamples
                        s_factor = 1/length(od.msn_away_lfr_dist);
                        choice = floor(rand()/s_factor)+1;
                        sub_idx = randsample(1:od.fsi_res.away_lfr_spec{iM}.spk_count, od.msn_away_lfr_dist(choice));
                        sub_idx = sort(sub_idx); sub_sts = this_sts;
                        sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                        sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                        sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                        temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                        sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                        temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                    end
                    od.fsi_res.away_lfr_spec{iM}.subsampled_sts = mean(temp_sts_vals,1);
                    od.fsi_res.away_lfr_spec{iM}.subsampled_ppc = mean(temp_ppc_vals,1);
                end
            end

             % Do control split stuff only if lfr/mfr splits are perfect
            % i.e., enough spikes and no nan flags

            if ~od.fsi_res.away_lfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.fsi_res.away_hfr_spec{iM}.flag_tooFewSpikes && ...
                    ~od.fsi_res.away_lfr_spec{iM}.flag_nansts && ...
                    ~od.fsi_res.away_lfr_spec{iM}.flag_nanppc  && ...
                    ~od.fsi_res.away_hfr_spec{iM}.flag_nansts && ...
                    ~od.fsi_res.away_hfr_spec{iM}.flag_nanppc && ...
                    od.fsi_res.away_spec{iM}.valid_split_count == cfg_master.nControlSplits

                od.fsi_res.away_spec{iM}.flag_no_control_split = false;

                % Create temporary arrays to average later
                all_sta_vals = zeros([2, cfg_master.nControlSplits, size(od.fsi_res.away_spec{iM}.sta_vals)]);
                all_freqs = zeros([cfg_master.nControlSplits, size(od.fsi_res.away_spec{iM}.freqs)]);
                all_sts_vals = zeros([2, cfg_master.nControlSplits, length(od.fsi_res.away_spec{iM}.sts_vals)]);
                all_ppc = zeros([2, cfg_master.nControlSplits, size(od.fsi_res.away_spec{iM}.ppc)]);
                all_spk_count = zeros(2, cfg_master.nControlSplits);
                all_subsampled_sts = zeros([2, cfg_master.nControlSplits, length(od.fsi_res.away_spec{iM}.sts_vals)]);
                all_subsampled_ppc = zeros([2, cfg_master.nControlSplits, size(od.fsi_res.away_spec{iM}.ppc)]);

                % For each split, fill in the tables
                % If any one of them is nan, exit loop, set
                % flag_no_control_split as true and move on to the next
                % cell    
                flag_nan_in_split = false;      
                for iSplit = 1:cfg_master.nControlSplits
                    % Calculate and Save all spec results for away_lfr_data
                    this_p1_trials = find(od.fsi_res.away_spec{iM}.valid_splits(iSplit,:));
                    this_p1_cfg.trl = cfg_away_trials.trl(this_p1_trials,:);
                    this_p1_data = ft_redefinetrial(this_p1_cfg, this_data);
                    this_p2_trials = find(~od.fsi_res.away_spec{iM}.valid_splits(iSplit,:));
                    this_p2_cfg.trl = cfg_away_trials.trl(this_p2_trials,:);
                    this_p2_data = ft_redefinetrial(this_p2_cfg, this_data);

                    p1_spk_count = 0;
                    for iT = 1:length(this_p1_data.trial)
                        p1_spk_count = p1_spk_count + sum(this_p1_data.trial{iT}(2,:)); 
                    end

                    p2_spk_count = 0;
                    for iT = 1:length(this_p2_data.trial)
                        p2_spk_count = p2_spk_count + sum(this_p2_data.trial{iT}(2,:)); 
                    end

                    all_spk_count(1,iSplit) = p1_spk_count;
                    all_spk_count(2,iSplit) = p2_spk_count;

                    % Calculate and save STA
                    this_p1_sta = ft_spiketriggeredaverage(cfg_ft, this_p1_data);
                    this_p2_sta = ft_spiketriggeredaverage(cfg_ft, this_p2_data);
                    all_sta_vals(1, iSplit, :) = this_p1_sta.avg(:,:)';
                    all_sta_vals(2, iSplit, :) = this_p2_sta.avg(:,:)';

                    % Calculate and save STS
                    cfg_sts.method = 'mtmconvol';
                    cfg_sts.foi = 1:1:100;
                    cfg_sts.t_ftimwin = 5./cfg_sts.foi;
                    cfg_sts.taper = 'hanning';
                    cfg_sts.spikechannel =  sd.S.ft_spikes(iC).label{1};
                    cfg_sts.channel = away_lfr_data.label{1};
                    this_p1_sts = ft_spiketriggeredspectrum(cfg_sts, this_p1_data);
                    this_p2_sts = ft_spiketriggeredspectrum(cfg_sts, this_p2_data);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_sts.fourierspctrm{1}),1)) || ...
                            ~isempty(find(isnan(this_p2_sts.fourierspctrm{1}),1))
                        flag_nan_in_split = true;
                        break;
                    end
                    all_freqs(iSplit,:) = this_p1_sts.freq;
                    all_sts_vals(1,iSplit,:) = nanmean(sq(abs(this_p1_sts.fourierspctrm{1})));
                    all_sts_vals(2,iSplit,:) = nanmean(sq(abs(this_p2_sts.fourierspctrm{1})));

                    % Calculate and save PPC
                    cfg_ppc               = [];
                    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg_ppc.spikechannel  = this_sts.label;
                    cfg_ppc.channel       = this_sts.lfplabel; % selected LFP channels
                    cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
                    cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
                    this_p1_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p1_sts);
                    this_p2_ppc              = ft_spiketriggeredspectrum_stat(cfg_ppc,this_p2_sts);
                    % Check for nans and abort if found
                    if ~isempty(find(isnan(this_p1_ppc.ppc0),1)) || ...
                         ~isempty(find(isnan(this_p2_ppc.ppc0),1))   
                        flag_nan_in_split = true;
                        break;
                    end
                    all_ppc(1,iSplit,:) = this_p1_ppc.ppc0';
                    all_ppc(2,iSplit,:) = this_p2_ppc.ppc0'; 

                    % Calculate and save subsampled measures for away Reward p1 trials
                    % When an FSI has fewer spikes than at least one MSN, do NOT
                    % subsample and include all the spikes
                    if isempty(od.msn_away_p1_dist)
						od.fsi_res.away_p1_spec{iM}.flag_no_subsampling = true;
                    elseif p1_spk_count < max(od.msn_away_p1_dist)
                        od.fsi_res.away_p1_spec{iM}.flag_no_subsampling = true;
                    elseif isfield(od.fsi_res.away_spec{iM},'flag_no_subsampling') && od.fsi_res.away_spec{iM}.flag_no_subsampling
                        od.fsi_res.away_p1_spec{iM}.flag_no_subsampling = true;
					else
                        od.fsi_res.away_p1_spec{iM}.flag_no_subsampling = false;
                            temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_spec{iM}.sts_vals));
                            temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_spec{iM}.ppc));
                        for iS = 1:cfg_master.num_subsamples
                            s_factor = 1/length(od.msn_away_p1_dist);
                            choice = floor(rand()/s_factor)+1;
                            sub_idx = randsample(1:p1_spk_count, od.msn_away_p1_dist(choice));
                            sub_idx = sort(sub_idx); 
                            sub_sts = this_p1_sts;
                            sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                            sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                            sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                            temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                            sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                            temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                        end
                        all_subsampled_sts(1,iSplit,:) = mean(temp_sts_vals,1);
                        all_subsampled_ppc(1,iSplit,:) = mean(temp_ppc_vals,1);
                    end

                    % Calculate and save subsampled measures for away Reward p2 trials
                    % When an FSI has fewer spikes than at least one MSN, do NOT
                    % subsample and include all the spikes
                    if isempty(od.msn_away_p2_dist)
						od.fsi_res.away_p2_spec{iM}.flag_no_subsampling = true;
                    elseif p2_spk_count < max(od.msn_away_p2_dist)
                        od.fsi_res.away_p2_spec{iM}.flag_no_subsampling = true;
                    elseif isfield(od.fsi_res.away_spec{iM},'flag_no_subsampling') && od.fsi_res.away_spec{iM}.flag_no_subsampling
                        od.fsi_res.away_p2_spec{iM}.flag_no_subsampling = true;
					else
                        od.fsi_res.away_p2_spec{iM}.flag_no_subsampling = false;
                            temp_sts_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_spec{iM}.sts_vals));
                            temp_ppc_vals = zeros(cfg_master.num_subsamples, length(od.fsi_res.away_spec{iM}.ppc));
                        for iS = 1:cfg_master.num_subsamples
                            s_factor = 1/length(od.msn_away_p2_dist);
                            choice = floor(rand()/s_factor)+1;
                            sub_idx = randsample(1:p2_spk_count, od.msn_away_p2_dist(choice));
                            sub_idx = sort(sub_idx); 
                            sub_sts = this_p2_sts;
                            sub_sts.fourierspctrm{1} = sub_sts.fourierspctrm{1}(sub_idx,:,:);
                            sub_sts.time{1} = sub_sts.time{1}(sub_idx,:);
                            sub_sts.trial{1} = sub_sts.trial{1}(sub_idx,:);
                            temp_sts_vals(iS,:) = nanmean(sq(abs(sub_sts.fourierspctrm{1})));
                            sub_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, sub_sts);
                            temp_ppc_vals(iS,:) = sub_ppc.ppc0';
                        end
                        all_subsampled_sts(2,iSplit,:) = mean(temp_sts_vals,1);
                        all_subsampled_ppc(2,iSplit,:) = mean(temp_ppc_vals,1);
                    end
                end
                if flag_nan_in_split
                   od.fsi_res.away_spec{iM}.flag_no_control_split = true;
                   break;
                end
                % Save averages and std
                od.fsi_res.away_p1_spec{iM}.sta_time = od.fsi_res.away_spec{iM}.sta_time;
                od.fsi_res.away_p2_spec{iM}.sta_time = od.fsi_res.away_spec{iM}.sta_time;
                od.fsi_res.away_p1_spec{iM}.freqs = squeeze(mean(all_freqs(1,:,:),2))';
                od.fsi_res.away_p2_spec{iM}.freqs = squeeze(mean(all_freqs(2,:,:),2))';
                od.fsi_res.away_p1_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(1,:,:),2))';
                od.fsi_res.away_p2_spec{iM}.mean_sta = squeeze(mean(all_sta_vals(2,:,:),2))';
                od.fsi_res.away_p1_spec{iM}.sd_sta = squeeze(std(all_sta_vals(1,:,:)))';
                od.fsi_res.away_p2_spec{iM}.sd_sta = squeeze(std(all_sta_vals(2,:,:)))';
                od.fsi_res.away_p1_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(1,:,:),2))';
                od.fsi_res.away_p2_spec{iM}.mean_sts = squeeze(mean(all_sts_vals(2,:,:),2))';
                od.fsi_res.away_p1_spec{iM}.sd_sts = squeeze(std(all_sts_vals(1,:,:)))';
                od.fsi_res.away_p2_spec{iM}.sd_sts = squeeze(std(all_sts_vals(2,:,:)))';
                od.fsi_res.away_p1_spec{iM}.mean_ppc = squeeze(mean(all_ppc(1,:,:),2))';
                od.fsi_res.away_p2_spec{iM}.mean_ppc = squeeze(mean(all_ppc(2,:,:),2))';
                od.fsi_res.away_p1_spec{iM}.sd_ppc = squeeze(std(all_ppc(1,:,:)))';
                od.fsi_res.away_p2_spec{iM}.sd_ppc = squeeze(std(all_ppc(2,:,:)))';
                od.fsi_res.away_p1_spec{iM}.mean_subsampled_sts = squeeze(mean(all_subsampled_sts(1,:,:),2))';
                od.fsi_res.away_p2_spec{iM}.mean_subsampled_sts = squeeze(mean(all_subsampled_sts(2,:,:),2))';
                od.fsi_res.away_p1_spec{iM}.sd_subsampled_sts = squeeze(std(all_subsampled_sts(1,:,:)))';
                od.fsi_res.away_p2_spec{iM}.sd_subsampled_sts = squeeze(std(all_subsampled_sts(2,:,:)))';
                od.fsi_res.away_p1_spec{iM}.mean_subsampled_ppc = squeeze(mean(all_subsampled_ppc(1,:,:),2))';
                od.fsi_res.away_p2_spec{iM}.mean_subsampled_ppc = squeeze(mean(all_subsampled_ppc(2,:,:),2))';
                od.fsi_res.away_p1_spec{iM}.sd_subsampled_ppc = squeeze(std(all_subsampled_ppc(1,:,:)))';
                od.fsi_res.away_p2_spec{iM}.sd_subsampled_ppc = squeeze(std(all_subsampled_ppc(2,:,:)))';
                od.fsi_res.away_p1_spec{iM}.spk_count = round(mean(all_spk_count(1,:)));
                od.fsi_res.away_p2_spec{iM}.spk_count = round(mean(all_spk_count(2,:)));
            else
                od.fsi_res.away_spec{iM}.flag_no_control_split = true;
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