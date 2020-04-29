%% temporary version to collect spike phases

clear
restoredefaultpath;
addpath(genpath('D:\My_Documents\GitHub\striatal-spike-rhythms\shared'));
addpath(genpath('D:\My_Documents\GitHub\striatal-spike-rhythms\chronux_2_12\spectral_analysis'));
%addpath(genpath('C:\Users\mvdm\Documents\GitHub\striatal-spike-rhythms\shared'));
%addpath(genpath('C:\Users\mvdm\Documents\GitHub\striatal-spike-rhythms\chronux_2_12\spectral_analysis'));

% set default interpreter to non-LaTeX

cfg_master = [];
cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
cfg_master.ccMethod = 'MvdM'; % method used for cell classification (see also 'Sirota')
cfg_master.minSpikes = 100; % only keep cells with at least this many spikes
cfg_master.plot = 0; % produce output figure for each cell?
cfg_master.wsize = 5; % length (s) of data analysis window either side of reward delivery
cfg_master.pad = 0.5; % padding (s) of data analysis window on either side (to avoid edge effects)
cfg_master.rats = {'R117', 'R119', 'R131', 'R132'};
cfg_master.nShuf = 0; % number of shuffles for spike-triggered spectrum and PPC
cfg_master.spk_dt = 0.0025; % interspike interval for surrogate spike train used for spike-triggered spectrum pool

[fd, fd_extra] = getDataPath(cfg_master);

%%
cc = 1; % cell count

%%
for iS = 1:length(fd)

    fprintf('Entering session %d/%d...\n',iS,length(fd));
    cd(fd{iS});
    
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
    
    % load reward events
    reward_t = getRewardTimes;
    fprintf('%d trials detected.\n', length(reward_t));
    
    if length(reward_t) < 20
        disp('*** SESSION SKIPPED DUE TO INSUFFICIENT TRIALS ***');
        continue;
    end
    
    %% load ft format data
    prev_path = set_ft_path;
    ft_lfp = ft_read_neuralynx_interp(ExpKeys.goodGamma_vStr);
    
    % awful time correction
    scale_fac = (csc.tvec(end)-csc.tvec(1)) ./ (ft_lfp.time{1}(end));
    ft_lfp.time{1} = ft_lfp.time{1} .* scale_fac;
    ft_lfp.hdr.Fs = round(ft_lfp.hdr.Fs ./ scale_fac);
    ft_lfp.fsample = ft_lfp.hdr.Fs;

    % figure out nlx-ft offset
    t0 = double(ft_lfp.hdr.FirstTimeStamp) ./ 10^6;

    reward_t = reward_t - t0;
    
    nan_idx = find(isnan(ft_lfp.trial{1}));
    ft_lfp.trial{1}(nan_idx) = 0;
    fprintf('  %d NaNs replaced in ft LFP.\n', length(nan_idx));
    
    % trialify data
    path(prev_path);
    reward_idx = nearest_idx3(reward_t, ft_lfp.time{1}); % now in samples
    prev_path = set_ft_path;
    
    dSamples = (cfg_master.wsize + cfg_master.pad) .* ft_lfp.hdr.Fs;
    start_idx = reward_idx - dSamples; end_idx = reward_idx + dSamples - 1;
    
    cfg_trl = []; cfg_trl.trl = cat(2, start_idx, end_idx, repmat(-dSamples, [length(start_idx) 1]));
    ft_lfp_trl = ft_redefinetrial(cfg_trl, ft_lfp);

    path(prev_path);
    
    %% load spike data
    cfg = []; cfg.uint = '32';
    S = LoadSpikes(cfg);
    S = restrict(S, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
    
    nSpikes = cellfun(@length, S.t);
    keep = nSpikes >= cfg_master.minSpikes; % note that there will be another exclusion step later (after restricting to analysis window)
    S = SelectTS([], S, keep);
    
    %% Categorize cells and add tetrode depths (would be nice to make this into a function)
    cfg_wv = []; cfg_wv.cMethod = cfg_master.ccMethod;
    s_out = CategorizeStriatumWave(cfg_wv, S);
    
    s_out.unit = [s_out.other s_out.msn s_out.fsi];
    s_out.ident = [zeros(1, length(s_out.other)) ones(1, length(s_out.msn)) repmat(2, 1, length(s_out.fsi))];
    
    cfg_tt = []; cfg_tt.verbose = 1;
    cfg_tt.this_rat_ID = fd_extra.ratID_num(iS);
    cfg_tt.this_date = fd_extra.fd_date_num(iS);
    
    for iC = 1:length(S.t)
         S.usr.cell_type(iC) = s_out.ident(find(s_out.unit == iC));
         S.usr.tetrodeDepths(iC) = ExpKeys.TetrodeDepths(S.usr.tt_num(iC));
         
         cfg_tt.ttno = S.usr.tt_num(iC);
         [S.usr.distanceTurned(iC), prev_fd] = DistanceTurned(cfg_tt, fd, fd_extra);
         cfg_tt.verbose = 0;
    end
    
    % correlate with previous session waveforms if available
    if isempty(prev_fd) % no previous day available
       S.usr.duplicate = zeros(size(S.usr.tt_num)); 
    else
       pushdir(prev_fd);
       S2 = LoadSpikes([]);
       nSpikes = cellfun(@length, S2.t); keep = nSpikes >= cfg_master.minSpikes;
       S2 = SelectTS([], S2, keep);
     
       s_out2 = CategorizeStriatumWave(cfg_wv, S2);
       s_out = CalcWVDistances([], s_out, s_out2); % add comparison with previous day's waveforms
       
       popdir;
                     
       % for each cell in current session, decide if duplicate
       for iC = 1:length(S.t)
          
           this_tt_no = S.usr.tt_num(iC);
           prev_day_cells = find(S2.usr.tt_num == this_tt_no);
           
           if isempty(prev_day_cells) % no cells recorded fron this tt in previous session
               S.usr.duplicate(iC) = 0;
           else % previous day cells found
               temp_corr = s_out.corr(iC, prev_day_cells);
               temp_peakn = s_out.peakdiffn(iC, prev_day_cells);
               
               if temp_corr > cfg_master.maxPrevCorr & abs(temp_peakn) < cfg_master.maxPeakn % wv correlation big, peak difference small
                   S.usr.duplicate(iC) = 1;
               else
                   S.usr.duplicate(iC) = 0;
               end
           end
       end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% main loop across cells %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iC = length(S.t):-1:1
        
        fprintf('Cell %d/%d...\n', iC, length(S.t));
        
        % skip if not turned & correlated with prev session
        if S.usr.duplicate(iC) & S.usr.distanceTurned(iC) < 80
            fprintf('Cell skipped - likely duplicate.\n');
            continue;
        end
        
        % skip if on same tt as LFP
        if S.usr.tt_num(iC) == lfp_tt
            fprintf('Cell skipped - same TT as LFP.\n');
            continue;
        end
        
        spk_t = S.t{iC};
        
        % make trials -- need to do this here because insufficient spikes trigger a skip
        clear tr_data spk_count;
        spk_count = 0;
        for iT = 1:length(reward_t) % this bit should be made into a function
            this_spk = spk_t - (reward_t(iT) - cfg_master.wsize - cfg_master.pad) - t0; % t0 correction necessary because reward_t is in ft timebase but spk_t is not :-(
            keep = this_spk >= 0 & this_spk < 2*cfg_master.wsize + 2*cfg_master.pad; spk_count = spk_count + sum(keep);
            tr_data(iT).times = this_spk(keep);
        end
        
        if spk_count < 100
            fprintf('Cell skipped - insufficient spikes.\n');
            continue;
        end
        
        % track some things about this cell
        ALL.cellType(cc) = S.usr.cell_type(iC);
        ALL.cellLabel{cc} = S.label(iC);
        ALL.cellDepth(cc) = S.usr.tetrodeDepths(iC);
        ALL.ID(cc) = s_out.cq.id(iC); ALL.Lr(cc) = s_out.cq.lr(iC); ALL.ampl(cc) = s_out.cq.ampl(iC);              
        ALL.sessno(cc) = iS;
        ALL.ratID(cc) = fd_extra.ratID_num(iS);
              
        %%%%%%%%%%%
        %%% ppc %%%
        %%%%%%%%%%%
        prev_path = set_ft_path;
        spike = ft_read_spike(S.label{iC}); % needs fixed read_mclust_t.m (edited D:\My_Documents\GitHub\fieldtrip\fileio\private\read_mclust_t.m to uint32)
                       
        % cut spike data in trials
        cfg           = [];
        cfg.hdr       = ft_lfp_trl.hdr; % contains information for conversion of samples to timestamps
        cfg.trlunit   = 'samples';
        cfg.trl       = cfg_trl.trl;
        spike_trl     = ft_spike_maketrials(cfg, spike);
        
        % spike triggered spectrum (nice convolution method)
        cfg              = [];
        cfg.method       = 'mtmconvol';
        cfg.foi          = 1:0.5:100;
        cfg.t_ftimwin    = 7./cfg.foi; % 5 cycles per frequency
        cfg.taper        = 'hanning';
        cfg.channel      = ft_lfp_trl.label(1);
        cfg.rejectsaturation = 'no';
        stsConvol        = ft_spiketriggeredspectrum(cfg, ft_lfp_trl, spike_trl);
        
        % compute the statistics on the phases
        cfg               = [];
        cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
        cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
        cfg.timwin        = 'all'; % compute over all available spikes in the window
        statSts           = ft_spiketriggeredspectrum_stat(cfg, stsConvol);
           
        ALL.ppc(cc,:) = statSts.ppc0';
        ALL.ppc_freq = statSts.freq;
        
        cfg.method        = 'ang';
        cfg.foi           = 'all';
        statSts           = ft_spiketriggeredspectrum_stat(cfg, stsConvol);
        
        ALL.ppc_ang(cc,:) = statSts.ang;
        
        cfg.method        = 'ral';
        cfg.foi           = 'all';
        statSts           = ft_spiketriggeredspectrum_stat(cfg, stsConvol);
        
        ALL.ppc_r(cc,:) = statSts.ral;
        
        if any(isnan(ALL.ppc(cc,:)))
           error('PPC failed!'); 
        end
        
        %%%%%%%%%%%
        %%% STA %%%
        %%%%%%%%%%%
        data_sta = ft_appendspike([], ft_lfp, spike_trl);
        
        cfg              = [];
        cfg.timwin       = [-2.5 2.5];
        cfg.spikechannel = spike_trl.label{1};
        cfg.channel      = ft_lfp.label{1};
        staAll           = ft_spiketriggeredaverage(cfg, data_sta);
               
        % add STA to ALL
        ALL.STA(cc,:) = staAll.avg; ALL.STAtime = staAll.time;
        
        % compute STA spectrum
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 0.5:0.5:100;
        cfg.t_ftimwin    = 7./cfg.foi;
        cfg.toi          = 0;
        staPow = ft_freqanalysis(cfg, staAll);
        
        ALL.STAp(cc,:) = squeeze(staPow.powspctrm); ALL.STAp_freq = staPow.freq;
        
        % compute STA phase
        cfg.output       = 'fourier';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 0.5:0.5:100;
        cfg.t_ftimwin    = 7./cfg.foi;
        cfg.toi          = 0;
        staPow = ft_freqanalysis(cfg, staAll);
        
        ALL.STAa(cc,:) = angle(squeeze(staPow.fourierspctrm));
        if any(ALL.STAa(cc,:) > pi)
            error('Impossible angle.')
        end
              
        path(prev_path);
 
        cc = cc + 1; % advance cell count
    end % of cell loop
    
end % of session loop


%%
function prev_path = set_ft_path()
% set path for FieldTrip

prev_path = path;

restoredefaultpath;
addpath('D:\My_Documents\GitHub\fieldtrip'); ft_defaults
addpath('D:\My_Documents\GitHub\striatal-spike-rhythms\shared\io\ft');
addpath('D:\My_Documents\GitHub\striatal-spike-rhythms\shared\io\neuralynx');

%addpath('C:\Users\mvdm\Documents\GitHub\fieldtrip'); ft_defaults
%addpath('C:\Users\mvdm\Documents\GitHub\striatal-spike-rhythms\shared\io\ft');
%addpath('C:\Users\mvdm\Documents\GitHub\striatal-spike-rhythms\shared\io\neuralynx');



end