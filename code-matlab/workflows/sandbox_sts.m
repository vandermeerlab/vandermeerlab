%% sandbox to generate trial-by-trial spike triggered spectrum of neurons
%
%
% The analysis proceeds as follows:
% - load data
% - restrict spikes to around +/-5 sec window around reward times and retain
%   neurons with at least 100 spikes remaining
% - separate type 1 (MSNs) and type 2 (FSI) neurons
% - generate trial-by-trial STS
%
% output variables are stored in the od variable, which can be saved and
% later used to plot the figures

%% setup
clear;
cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/');
please = [];
please.rats = {'R117'}%,'R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.output_dir = '/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp';
cfg_in.cx_binsize = 0.01; 

%%
% Top level loop which calls the main function for all the sessions
for iS = 1%:length(cfg_in.fd) % for each session...
    
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
    cfg_master.output_prefix = 'ccf_';
    cfg_master.output_dir = 'C:\temp';
    cfg_master.exc_types = 0; %cell types to be excluded
    cfg_master.cx_msn = 0.01;  % bin size for cross-correlation between pairs of Medium Spiny Neurons (MSNs) in seconds
    cfg_master.cx_fsi = 0.005;  % bin size for cross-correlations between pairs of Fast Spiking Interneurons (FSIs) in seconds
    cfg_master.cx_mix = 0.005;  % bin size for cross-correlations between pairs of 1 FSI abd 1 MSN in seconds
    cfg_master.max_t = 0.5; % half window length for the cross correlations in seconds in seconds
    cfg_master.plotfft = 1;
    cfg_master.plot = 1;

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
    % times
    rt = getRewardTimes();
    trial_starts = rt - 5;
    trial_ends = rt + 5;
    rt_iv = iv(trial_starts, trial_ends);
    rt_iv = MergeIV([],rt_iv);
    od.S = restrict(sd.S, rt_iv);
    
    od.S.cell_type = sd.S.usr.cell_type;
    od.S.tt_id = sd.S.usr.tt_num;

    % Keep cells greater with greater than nMinSpike spikes and of the allowed types
    od.S = KeepCells(od.S,cfg_master.nMinSpikes,cfg_master.exc_types,lfp_tt);
    
    % parameters for STS
    cfg_s = []; cfg_s.binsize = 0.001; cfg_s.Fs = 1 ./ cfg_s.binsize;
%     cfg_s.nShuf = cfg_master.nShuf; 
    cfg_s.trial_len = 10;
    
    cfg.params = []; cfg.params.Fs = 200; cfg.params.tapers = [3 5];
    

    for iC = length(od.S.t):-1:1
         S = SelectTS([],od.S,iC);
        %calculating trial-trial STA
        for iT = 1:length(trial_starts)
            triaL_iv = iv(trial_starts(iT),trial_ends(iT));
            trial_S = restrict(S, triaL_iv);
            spk_t = trial_S.t{1};
            % Bin spikes on CSC timebase
            this_spk_binned = zeros(size(csc.tvec));
            idx = nearest_idx3(spk_t, csc.tvec);
            this_spk_binned(idx) = 1;
            [xc, tvec] = xcorr(csc.data, this_spk_binned, 5000);
            tvec = tvec .* median(diff(csc.tvec));

            if cfg_master.plotfft
                subplot(211);
                plot(tvec, xc);
                title('raw STA');
            end
% 
            % STA spectrum
            this_Fs = 1 ./ median(diff(csc.tvec));
            [P,F] = pwelch(xc,length(xc), [], [], this_Fs);

            if cfg_master.plot
                subplot(212);
                plot(F, 10*log10(P)); xlim([0 100]);
                title('STA spectrum');
            end
        end
    end
end

% Other functions
%%
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
