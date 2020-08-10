%% script to generate STA spectra and average STS on a trial-by trial basis as well as for binned trials
% trials are binned on the basis of mean firing rate in the trial such that
% the number of spikes in each of the bins are as close as possible
% The spectra for the FSIs are calculated after multiple rounds of
% subsampling
%% setup
clear;
cd('D:\ADRLabData');
please = [];
please.rats = {'R117','R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.output_dir = 'D:\RandomVstrAnalysis\temp';
cfg_in.exc_types = 0;


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
    
    lfp_tt = regexp(cfg.fc, 'CSC\d+', 'match');
    lfp_tt = str2double(lfp_tt{1}{1}(4:end)); % need this to skip cells from same tt (could make into function)
    fprintf('LFP ttno is %d\n', lfp_tt);
    
    % params
    cfg_master = []; % overall params
    cfg_master.dt = 0.001;
    cfg_master.ccMethod = 'MvdM'; % cell type classification method
    cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
    cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
    cfg_master.nMinSpikes = 100;
    cfg_master.iS = 1; % current session number out of fd list, get this from input cfg
    cfg_master.fd = []; % full list of session fd's, get this from input cfg
    cfg_master.fd_extra = []; % get this from input cfg
    cfg_master.write_output = 0;
    cfg_master.output_prefix = 'sts_';
    cfg_master.output_dir = 'C:\temp';
    cfg_master.exc_types = 0; %cell types to be excluded

    cfg_master = ProcessConfig(cfg_master,cfg_in);
    
    % spikes
    sd.S = LoadSpikesTarget(cfg_master);

    % Categorize cells and add tetrode depths
    cfg_wv = []; cfg_wv.cMethod = cfg_master.ccMethod;
    s_out = CategorizeStriatumWave(cfg_wv, sd.S);

    s_out.unit = [s_out.other s_out.msn s_out.fsi];
    s_out.ident = [zeros(1, length(s_out.other)) ones(1, length(s_out.msn)) repmat(2, 1, length(s_out.fsi))];

    cfg_tt = []; cfg_tt.verbose = 1;
    cfg_tt.this_rat_ID = cfg_master.fd_extra.ratID_num(cfg_master.iS);
    cfg_tt.this_date = cfg_master.fd_extra.fd_date_num(cfg_master.iS);

    for iC = 1:length(sd.S.t)
        sd.S.usr.cell_type(iC) = s_out.ident(find(s_out.unit == iC));
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
        nSpikes = cellfun(@length, S2.t); keep = nSpikes >= cfg_master.nMinSpikes;
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
    rt_iv = MergeIV([], rt_iv);
    % Saving near trial fields
    od.S1 = restrict(sd.S, rt_iv);
    od.S1.cell_type = sd.S.usr.cell_type;
    od.S1.tt_id = sd.S.usr.tt_num;
    od.S1.trial_starts = w_start;
    od.S1.trial_ends = w_end;
    
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
    rt_iv = MergeIV([], rt_iv);
    % Saving away trial fields
    od.S2 = restrict(sd.S, rt_iv);
    od.S2.cell_type = sd.S.usr.cell_type;
    od.S2.tt_id = sd.S.usr.tt_num;
    od.S2.trial_starts = w_start;
    od.S2.trial_ends = w_end;

    % Keep cells greater with greater than nMinSpike spikes and of the allowed types
    od.S1 = KeepCells(od.S1,cfg_master.nMinSpikes,cfg_master.exc_types,lfp_tt);
    od.S2 = KeepCells(od.S2,cfg_master.nMinSpikes,cfg_master.exc_types,lfp_tt);
    
    % Parameters for STS
    % Use dummy data to get Frequency values of pwelch and mtspectumFc
    cfg_s.Fs = 1/median(diff(csc.tvec));
    cfg_s.sts_wl = round(cfg_s.Fs);
    dd = zeros(1,cfg_s.sts_wl);
    [~,F] = pwelch(dd,length(dd),[],[],cfg_s.Fs); 
    cfg_s.foi = (F >= 0 & F <= 100);
    cfg_s.f_len = sum(cfg_s.foi);
    cfg_s.cfg_mt.tapers = [3 5];
    cfg_s.cfg_mt.Fs = cfg_s.Fs;
    [~,F] = mtspectrumc(dd,cfg_s.cfg_mt);
    cfg_s.foi_mt = (F >= 0 & F <= 100);
    cfg_s.f_len_mt = sum(cfg_s.foi_mt);
    cfg_s.lfp_data = csc.data;
    cfg_s.lfp_ts = csc.tvec;
    cfg_s.pbins = 2;
    cfg_s.num_samples = 1000;
    
    % determining spike triggered LFP segment length using dummy data
    [dum_seg,~]  = xcorr(csc.data, csc.data, floor(cfg_s.sts_wl/2));
    cfg_s.seg_len = length(dum_seg);
    tcount = length(od.S1.trial_starts);
    od.S1.freqs = F(F >= 0 & F <= 100);
    od.S1.msn_res = repmat(struct('sta_mtspec_ptile', zeros(cfg_s.pbins, cfg_s.f_len_mt), ...
                    'sta_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'mtsts_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'scount_ptile',zeros(cfg_s.pbins, 1), ...
                    'mfr', zeros(1, tcount), ...
                    'ptile_mfrs', zeros(cfg_s.pbins,2)), length(find(od.S1.cell_type == 1)), 1);
    od.S1.fsi_res = repmat(struct('sta_mtspec_ptile', zeros(cfg_s.pbins, cfg_s.f_len_mt), ...
                    'sta_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'mtsts_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'scount_ptile',zeros(cfg_s.pbins, 1), ...
                    'scount_sub_ptile',zeros(cfg_s.pbins, cfg_s.num_samples), ...
                    'mfr', zeros(1, tcount), ...
                    'ptile_mfrs', zeros(cfg_s.pbins,2)), length(find(od.S1.cell_type == 2)), 1);
    tcount = length(od.S2.trial_starts);           
    od.S2.freqs = F(F >= 0 & F <= 100);
    od.S2.msn_res = repmat(struct('sta_mtspec_ptile', zeros(cfg_s.pbins, cfg_s.f_len_mt), ...
                    'sta_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'mtsts_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'scount_ptile',zeros(cfg_s.pbins, 1), ...
                    'mfr', zeros(1, tcount), ...
                    'ptile_mfrs', zeros(cfg_s.pbins,2)), length(find(od.S2.cell_type == 1)), 1);
    od.S2.fsi_res = repmat(struct('sta_mtspec_ptile', zeros(cfg_s.pbins, cfg_s.f_len_mt), ...
                    'sta_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'mtsts_ptile', zeros(cfg_s.pbins, cfg_s.f_len), ...
                    'scount_ptile',zeros(cfg_s.pbins, 1), ...
                    'scount_sub_ptile',zeros(cfg_s.pbins, cfg_s.num_samples), ...
                    'mfr', zeros(1, tcount), ...
                    'ptile_mfrs', zeros(cfg_s.pbins,2)), length(find(od.S2.cell_type == 2)), 1);
                
    % For near trials
    msn = find(od.S1.cell_type == 1);
    msn_lfr_dist = zeros(length(msn),1);
    msn_hfr_dist = zeros(length(msn),1);
    for iC = 1:length(msn)
        S1 = SelectTS([],od.S1,msn(iC));
        cfg_s.trial_starts = od.S1.trial_starts;
        cfg_s.trial_ends = od.S1.trial_ends;
        od.S1.msn_res(iC) = calculateMSNspec(cfg_s, S1);
        msn_lfr_dist(iC) = od.S1.msn_res(iC).scount_ptile(1);
        msn_hfr_dist(iC) = od.S1.msn_res(iC).scount_ptile(2);
    end
    keep = (msn_lfr_dist >= 200) & (msn_hfr_dist >= 200);
    msn_lfr_dist = msn_lfr_dist(keep);
    msn_hfr_dist = msn_hfr_dist(keep);
    if numel(msn_lfr_dist) > 0
        fsi = find(od.S1.cell_type == 2);
        for iC = 1:length(fsi)
            S1 = SelectTS([],od.S1,fsi(iC));
            cfg_s.trial_starts = od.S1.trial_starts;
            cfg_s.trial_ends = od.S1.trial_ends;
            cfg_s.ldist = msn_lfr_dist;
            cfg_s.hdist = msn_hfr_dist;
            od.S1.fsi_res(iC) = calculateFSIspec(cfg_s, S1);
        end
    end
    % For away trials
    msn = find(od.S2.cell_type == 1);
    msn_lfr_dist = zeros(length(msn),1);
    msn_hfr_dist = zeros(length(msn),1);
    for iC = 1:length(msn)
        S2 = SelectTS([],od.S2,msn(iC));
        cfg_s.trial_starts = od.S2.trial_starts;
        cfg_s.trial_ends = od.S2.trial_ends;
        od.S2.msn_res(iC) = calculateMSNspec(cfg_s, S2);
        msn_lfr_dist(iC) = od.S2.msn_res(iC).scount_ptile(1);
        msn_hfr_dist(iC) = od.S2.msn_res(iC).scount_ptile(2);
    end
    keep = (msn_lfr_dist >= 200) & (msn_hfr_dist >= 200);
    msn_lfr_dist = msn_lfr_dist(keep);
    msn_hfr_dist = msn_hfr_dist(keep);
    if numel(msn_lfr_dist) > 0
        fsi = find(od.S2.cell_type == 2);
        for iC = 1:length(fsi)
            S2 = SelectTS([],od.S2,fsi(iC));
            cfg_s.trial_starts = od.S2.trial_starts;
            cfg_s.trial_ends = od.S2.trial_ends;
            cfg_s.ldist = msn_lfr_dist;
            cfg_s.hdist = msn_hfr_dist;
            od.S2.fsi_res(iC) = calculateFSIspec(cfg_s, S2);
        end
    end
    
    if cfg_master.write_output
         [~, fp, ~] = fileparts(pwd);
         pushdir(cfg_master.output_dir);
         fn_out = cat(2, fp, '_sts.mat');
         save(fn_out,'od'); % should add option to save in specified output dir
         popdir;
    end
    
end

%% Other functions

%% function that calculates STA spectrum  as well as the average spectra 
% for spike-triggered LFP using mutli-taper method for MSNs
function res = calculateMSNspec(cfg_in, S)
        tcount = length(cfg_in.trial_starts);
        res.mfr = zeros(tcount,1);
        all_tspikes = cell(1,tcount);
        for iT = 1:tcount
            triaL_iv = iv(cfg_in.trial_starts(iT),cfg_in.trial_ends(iT));
            trial_S = restrict(S, triaL_iv);
            all_tspikes{iT} = trial_S.t{1};
            res.mfr(iT) = length(all_tspikes{iT})/(cfg_in.trial_ends(iT) - cfg_in.trial_starts(iT));
        end
        
        % Get rid of trials with no spikes and bin the rest
        nz_trials = find(res.mfr ~= 0);
        nz_mfr = res.mfr(nz_trials);
        nz_tcount = length(nz_trials);
        nz_tspikes = cell(nz_tcount,1);
        spk_tcount = zeros(nz_tcount,1);
        for iT = 1:nz_tcount
            nz_tspikes{iT} = all_tspikes{nz_trials(iT)};
            spk_tcount(iT) = length(nz_tspikes{iT});
        end
        
        % Find out the firing rate threhsold to split the trials such that
        % spikes are more or less equally divided
        ufr = unique(nz_mfr);
        dif_min = sum(cell2mat(nz_tspikes));

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
        
        % Bin trials based on the firing rate threshold 
        % and calculate results for the bins
        res.ptile_mfrs = zeros(cfg_in.pbins,2);
        res.sta_ptile = zeros(cfg_in.pbins, cfg_in.sts_wl);
        res.sta_mtspec_ptile = zeros(cfg_in.pbins, cfg_in.f_len_mt);
        res.mtsts_ptile = zeros(cfg_in.pbins, cfg_in.f_len_mt);
        res.scount_ptile = zeros(cfg_in.pbins, 1);
        for iP = 1:cfg_in.pbins
            if (iP == 1)
                res.ptile_mfrs(iP,1) = min(nz_mfr);
                res.ptile_mfrs(iP,2) = fr_thresh;
                ptrials = find(nz_mfr >= res.ptile_mfrs(iP,1) & nz_mfr <= res.ptile_mfrs(iP,2));
            else
                res.ptile_mfrs(iP,1) = fr_thresh;
                res.ptile_mfrs(iP,2) = max(nz_mfr);
                ptrials = find(nz_mfr > res.ptile_mfrs(iP,1) & nz_mfr <= res.ptile_mfrs(iP,2));
            end
            p_spikes = cell(length(ptrials),1);
            if ~isempty(ptrials)
                for iT = 1:length(ptrials)
                    p_spikes{iT} = nz_tspikes{ptrials(iT)};     
                end
                p_spikes = cell2mat(p_spikes);
            end
            res.scount_ptile(iP) = length(p_spikes);
            % Calculate average STS of the mfr percentile based binned spikes
            % Bin spikes on CSC timebase
            idx = nearest_idx3(p_spikes, cfg_in.lfp_ts);
            this_p2 = zeros(length(idx),cfg_in.f_len_mt);
            % For each spike
            for iS = 1:length(idx)
                this_seg = cfg_in.lfp_data((idx(iS)-floor(cfg_in.sts_wl/2)):(idx(iS)+floor(cfg_in.sts_wl/2)));
                [P2,~] = mtspectrumc(this_seg, cfg_in.cfg_mt);
                this_p2(iS,:) = P2(cfg_in.foi_mt);
            end
            res.mtsts_ptile(iP,:) = mean(this_p2,1);
            % Calculate spectrum of STA for the mfr percentile based binned spikes
            this_spk_binned = zeros(size(cfg_in.lfp_ts));
            this_spk_binned(idx) = 1;
            [this_sta, ~] = xcorr(cfg_in.lfp_data, this_spk_binned, floor(cfg_in.sts_wl/2));
            [P2,~] = mtspectrumc(this_sta, cfg_in.cfg_mt);
            res.sta_mtspec_ptile(iP,:) = P2(cfg_in.foi_mt);
            res.sta_ptile(iP,:) = this_sta;
        end 
end


%% function that calculates STA spectrum  as well as the average spectra 
% for spike-triggered LFP using mutli-taper method for FSIs
function res = calculateFSIspec(cfg_in, S)
        tcount = length(cfg_in.trial_starts);
        res.mfr = zeros(tcount,1);
        all_tspikes = cell(1,tcount);
        for iT = 1:tcount
            triaL_iv = iv(cfg_in.trial_starts(iT),cfg_in.trial_ends(iT));
            trial_S = restrict(S, triaL_iv);
            all_tspikes{iT} = trial_S.t{1};
            res.mfr(iT) = length(all_tspikes{iT})/(cfg_in.trial_ends(iT) - cfg_in.trial_starts(iT));
        end
        
        % Get rid of trials with no spikes and bin the rest
        nz_trials = find(res.mfr ~= 0);
        nz_mfr = res.mfr(nz_trials);
        nz_tcount = length(nz_trials);
        nz_tspikes = cell(nz_tcount,1);
        spk_tcount = zeros(nz_tcount,1);
        for iT = 1:nz_tcount
            nz_tspikes{iT} = all_tspikes{nz_trials(iT)};
            spk_tcount(iT) = length(nz_tspikes{iT});
        end
        
        % Find out the firing rate threhsold to split the trials such that
        % spikes are more or less equally divided
        ufr = unique(nz_mfr);
        dif_min = sum(cell2mat(nz_tspikes));

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
        
        % Bin trials based on the firing rate threshold 
        % and calculate results for the bins
        res.ptile_mfrs = zeros(cfg_in.pbins,2);
        res.sta_ptile = zeros(cfg_in.pbins, cfg_in.sts_wl);
        res.sta_mtspec_ptile = zeros(cfg_in.pbins, cfg_in.f_len_mt);
        res.mtsts_ptile = zeros(cfg_in.pbins, cfg_in.f_len_mt);
        res.scount_ptile = zeros(cfg_in.pbins, 1);
        res.scount_sub_ptile = zeros(cfg_in.pbins, cfg_in.num_samples);
        for iP = 1:cfg_in.pbins
            if (iP == 1)
                res.ptile_mfrs(iP,1) = min(nz_mfr);
                res.ptile_mfrs(iP,2) = fr_thresh;
                ptrials = find(nz_mfr >= res.ptile_mfrs(iP,1) & nz_mfr <= res.ptile_mfrs(iP,2));
            else
                res.ptile_mfrs(iP,1) = fr_thresh;
                res.ptile_mfrs(iP,2) = max(nz_mfr);
                ptrials = find(nz_mfr > res.ptile_mfrs(iP,1) & nz_mfr <= res.ptile_mfrs(iP,2));
            end
            p_spikes = cell(length(ptrials),1);
            if ~isempty(ptrials)
                for iT = 1:length(ptrials)
                    p_spikes{iT} = nz_tspikes{ptrials(iT)};     
                end
                p_spikes = cell2mat(p_spikes);
            end
            res.scount_ptile(iP) = length(p_spikes);
            
            % When an FSI has fewer spikes than at least one MSN, do NOT
            % subsample and include all the spikes
            if length(p_spikes) < max(cfg_in.ldist) || length(p_spikes) < max(cfg_in.hdist)
                idx = nearest_idx3(p_spikes, cfg_in.lfp_ts);
                this_p2 = zeros(length(idx),cfg_in.f_len_mt);
                % For each spike
                for iS = 1:length(idx)
                    this_seg = cfg_in.lfp_data((idx(iS)-floor(cfg_in.sts_wl/2)):(idx(iS)+floor(cfg_in.sts_wl/2)));
                    [P2,~] = mtspectrumc(this_seg, cfg_in.cfg_mt);
                    this_p2(iS,:) = P2(cfg_in.foi_mt);
                end
                res.mtsts_ptile(iP,:) = mean(this_p2,1);
                % Calculate spectrum of STA for the mfr percentile based binned spikes
                this_spk_binned = zeros(size(cfg_in.lfp_ts));
                this_spk_binned(idx) = 1;
                [this_sta, ~] = xcorr(cfg_in.lfp_data, this_spk_binned, floor(cfg_in.sts_wl/2));
                [P2,~] = mtspectrumc(this_sta, cfg_in.cfg_mt);
                res.sta_mtspec_ptile(iP,:) = P2(cfg_in.foi_mt);
                res.sta_ptile(iP,:) = this_sta;
                continue
            end
            
            % Subsample spikes based on the MSN distribution and calculate
            % all the metrics for each sampling
            % TODO:
            s_factor = 1/length(cfg_in.ldist);
            as_sta = zeros(cfg_in.num_samples, cfg_in.sts_wl);
            as_mtspec = zeros(cfg_in.num_samples, cfg_in.f_len_mt);
            as_mtsts = zeros(cfg_in.num_samples, cfg_in.f_len_mt);
            for iR = 1:cfg_in.num_samples
                disp(strcat(2,'Subsampling round: ',num2str(iR)));
                % subsampling step
                choice = floor(rand()/s_factor)+1;
                if (iP == 1)
                    res.scount_sub_ptile(iP,iR) = cfg_in.ldist(choice);
                else
                    res.scount_sub_ptile(iP,iR) = cfg_in.hdist(choice);
                end
                this_spikes = randsample(p_spikes, res.scount_sub_ptile(iP,iR));
                
                idx = nearest_idx3(this_spikes, cfg_in.lfp_ts);
                this_p2 = zeros(length(idx),cfg_in.f_len_mt);
                % For each spike
                for iS = 1:length(idx)
                    this_seg = cfg_in.lfp_data((idx(iS)-floor(cfg_in.sts_wl/2)):(idx(iS)+floor(cfg_in.sts_wl/2)));
                    [P2,~] = mtspectrumc(this_seg, cfg_in.cfg_mt);
                    this_p2(iS,:) = P2(cfg_in.foi_mt);
                end
                as_mtsts(iR,:) = mean(this_p2,1);
                % Calculate spectrum of STA for the mfr percentile based binned spikes
                this_spk_binned = zeros(size(cfg_in.lfp_ts));
                this_spk_binned(idx) = 1;
                [this_sta, ~] = xcorr(cfg_in.lfp_data, this_spk_binned, floor(cfg_in.sts_wl/2));
                [P2,~] = mtspectrumc(this_sta, cfg_in.cfg_mt);
                as_mtspec(iR,:) = P2(cfg_in.foi_mt);
                as_sta(iR,:) = this_sta;
            end
            %Saving the mean of all the sampling rounds
            res.mtsts_ptile(iP,:) = mean(as_mtsts,1);
            res.sta_mtspec_ptile(iP,:) = mean(as_mtspec,1);
            res.sta_ptile(iP,:) = mean(as_sta,1);
        end 
end

%% function to remove cells based on spike threshold, cell type and tetrode
function S = KeepCells(S,minSpikes,cTypes,tt_num)
    thr = minSpikes;

    nCells = length(S.t);
    for i = nCells:-1:1
       l(i) = length(S.t{i}); 
    end
    keep = l >= thr;
    for iC = 1:length(keep)
        for iT = 1:length(cTypes)
            if S.cell_type(iC) == cTypes(iT)
                keep(iC) = 0;
                break
            end
        end
    end
    keep2 = (S.tt_id ~= tt_num);
    keep = keep & keep2;
    S.label = S.label(keep);
    S.t = S.t(keep);
    S.cell_type = S.cell_type(keep);
    S.tt_id = S.tt_id(keep);
end


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

