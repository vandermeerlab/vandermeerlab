%% sandbox to generate trial-by-trial average spike triggered spectrum of neurons
%
%
% The analysis proceeds as follows:
% - load data
% - separate type 1 (MSNs) and type 2 (FSI) neurons
% - spearate spikes into two epochs: near reward times and away from it
% - generate trial-by-trial STS
%
% output variables are stored in the od variable, which can be saved and
% later used to plot the figures

%% setup
clear;
cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/');
please = [];
please.rats = {'R117','R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.output_dir = '/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp';
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

    % Divide Spikes into two kinds of epochs: Around +/-5 s of the reward 
    % times (near trials) and the others (away trials).
    rt = getRewardTimes();
    near_trial_starts = rt - 5;
    near_trial_ends = rt + 5;
    rt_iv = iv(near_trial_starts, near_trial_ends);
    rt_iv = MergeIV([],rt_iv);
    od.S1 = restrict(sd.S, rt_iv);
    
    rt = getRewardTimes();
    away_trial_starts = [ExpKeys.TimeOnTrack;rt(1:end-1)+5];
    away_trial_ends = (rt-5);
    % Get rid of extra long-trials (greater than mean + 1*SD)
    trial_length = away_trial_ends - away_trial_starts;
    mtl = mean(trial_length); stl = std(trial_length);
    valid_trials = (trial_length <= mtl + stl); 
    away_trial_starts = away_trial_starts(valid_trials);
    away_trial_ends = away_trial_ends(valid_trials);
    rt_iv = iv(away_trial_starts, away_trial_ends);
    rt_iv = MergeIV([],rt_iv);
    od.S2 = restrict(sd.S, rt_iv);
    
    % Saving other specific fields
    od.S1.cell_type = sd.S.usr.cell_type;
    od.S1.tt_id = sd.S.usr.tt_num;
    od.S1.trial_starts = near_trial_starts;
    od.S1.trial_ends = near_trial_ends;
    od.S2.cell_type = sd.S.usr.cell_type;
    od.S2.tt_id = sd.S.usr.tt_num;
    od.S2.trial_starts = away_trial_starts;
    od.S2.trial_ends = away_trial_ends;

    % Keep cells greater with greater than nMinSpike spikes and of the allowed types
    od.S1 = KeepCells(od.S1,cfg_master.nMinSpikes,cfg_master.exc_types,lfp_tt);
    od.S2 = KeepCells(od.S2,cfg_master.nMinSpikes,cfg_master.exc_types,lfp_tt);
    
    % Parameters for STS
    % Use dummy data to get Frequency values of pwelch
    cfg_s.Fs = 1/median(diff(csc.tvec));
    cfg_s.sts_wl = round(cfg_s.Fs);
    dd = zeros(1,cfg_s.sts_wl);
    [~,F] = pwelch(dd,length(dd),[],[],cfg_s.Fs);
    cfg_s.foi = (F >= 0 & F <= 100);
    cfg_s.f_len = sum(cfg_s.foi);
    cfg_s.lfp_data = csc.data;
    cfg_s.lfp_ts = csc.tvec;
    
    od.S1.sts_p = cell(length(od.S1.t),1);
    od.S1.mfr = cell(length(od.S1.t),1);
    od.S2.sts_p = cell(length(od.S2.t),1);
    od.S2.mfr = cell(length(od.S2.t),1);
    
    % For near trials
    for iC = length(od.S1.t):-1:1 
        S1 = SelectTS([],od.S1,iC);
        [od.S1.sts{iC}, od.S1.mfr{iC}] = calculateSTS(cfg_s, S1, od.S1.trial_starts, od.S1.trial_ends);
    end
    
    % For away trials
    for iC = length(od.S2.t):-1:1 
        S2 = SelectTS([],od.S2,iC);
        [od.S2.sts{iC}, od.S2.mfr{iC}] = calculateSTS(cfg_s, S2, od.S2.trial_starts, od.S2.trial_ends);
    end
    
    if cfg_master.write_output
         [~, fp, ~] = fileparts(pwd);
         pushdir(cfg_master.output_dir);
         fn_out = cat(2,cfg_master.output_prefix, fp, '_sts.mat');
         save(fn_out,'od'); % should add option to save in specified output dir
         popdir;
    end
    
end

%% Other functions

%% function that calculates average STS for each trial
function [this_sts, this_mfr] = calculateSTS(cfg_in, S, trial_starts, trial_ends)
    this_Fs = cfg_in.Fs;
    sts_wl = cfg_in.sts_wl;
    f_len = cfg_in.f_len;
    lfp_data = cfg_in.lfp_data;
    lfp_ts = cfg_in.lfp_ts;
    foi = cfg_in.foi;
    % Calculating average spectra of individual STAs
    this_sts = zeros(length(trial_starts),f_len);
    this_mfr = zeros(length(trial_starts),1);
%     close all;
%     figure;
    for iT = 1:length(trial_starts)
        triaL_iv = iv(trial_starts(iT),trial_ends(iT));
        trial_S = restrict(S, triaL_iv);
        spk_t = trial_S.t{1};
        this_mfr(iT) = length(spk_t)/(trial_ends(iT) - trial_starts(iT));
        % Bin spikes on CSC timebase
        idx = nearest_idx3(spk_t, lfp_ts);
        this_p = zeros(length(idx),f_len);
        % For each spike
        for iS = 1:length(idx)
            % Using the fact that these spikes are well into the middle of
            % the session, ideally a function to pad with zeroes for the
            % STA if need be
            this_STA = lfp_data((idx(iS)-floor(sts_wl/2)):(idx(iS)+floor(sts_wl/2)));
            [P,~] = pwelch(this_STA,length(this_STA), [], [], this_Fs);
            this_p(iS,:) = P(foi);
        end
        this_sts(iT,:) = mean(this_p,1);
%         plot(10*log(this_sts(iT,:)), 'Color','blue');
%         hold on;
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
