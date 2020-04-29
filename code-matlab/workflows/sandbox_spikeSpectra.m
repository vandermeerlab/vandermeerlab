%% script exploring spike autocorrelations, spike spectra, and spike-triggered averages
% TODO:
% - add SFC and/or PPC
% - use alternate LFP if tt is same as LFP

% remember to set path
clear

% set default interpreter to non-LaTeX

cfg_master = [];
cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
cfg_master.maxPeakn = 0.1; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
cfg_master.ccMethod = 'MvdM'; % method used for cell classification (see also 'Sirota')
cfg_master.minSpikes = 100; % only keep cells with at least this many spikes
cfg_master.plot = 0; % produce output figure for each cell?
cfg_master.nShuf = 25; % number of spike shuffles
cfg_master.trial_len = 50; % length (s) of trials to chop data up into
cfg_master.rats = {'R117', 'R119', 'R131', 'R132'};

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

    %% load data
    cfg = []; cfg.uint = '32';
    S = LoadSpikes(cfg);
    S = restrict(S, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
    
    nSpikes = cellfun(@length, S.t);
    keep = nSpikes >= cfg_master.minSpikes;
    S = SelectTS([], S, keep);
    
    %% Categorize cells and add tetrode depths
    cfg_wv = [];cfg_wv.cMethod = cfg_master.ccMethod;
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

    %% Updated version with spike spectrum and STA
    % parameters for spike spectrum
    cfg = []; cfg.binsize = 0.001; cfg.Fs = 1 ./ cfg.binsize;
    cfg.nShuf = cfg_master.nShuf; cfg.trial_len = cfg_master.trial_len;
    
    cfg.params = []; cfg.params.Fs = 200; cfg.params.tapers = [3 5];
    
    tbin_edges = firstSpike(S):cfg.binsize:lastSpike(S);
    trial_starts = tbin_edges(1):cfg.trial_len:tbin_edges(end);
    
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
        
        % binarize spike train for acorr (display only)
        spk_binned = histc(spk_t, tbin_edges); spk_binned = spk_binned(1:end-1);
        [acf, tvec] = ComputeACF(cfg, spk_binned);
        
        if cfg_master.plot
            figure; subplot(221);
            plot(tvec,acf); xlim([-0.5 0.5]);
            title(sprintf('Cell %s', ALL_cellLabel{cc}));
        end
        
        % spike spectrum
        clear data;
        for iT = 1:length(trial_starts) % this bit should be made into a function
            this_spk = spk_t - trial_starts(iT);
            keep = this_spk >= 0 & this_spk < cfg.trial_len;
            data(iT).times = this_spk(keep);
        end
        
        [P,F,R] = mtspectrumpt(data, cfg.params);
        P = nanmean(P, 2) ./ length(spk_t);
        keep = F > 1;
        P = 10*log10(P(keep)); F = F(keep);
        
        if cfg_master.plot
            subplot(223);
            plot(F, medfilt1(P, 201), 'k', 'LineWidth', 2);
            drawnow;
        end
        
        if length(F) ~= 8111
            fprintf('Cell skipped - size mismatch.\n'); % unclear why Chronux sometimes returns different size
            continue;
        end
        
        % track some things about this cell
        ALL.cellType(cc) = S.usr.cell_type(iC);
        ALL.cellLabel{cc} = S.label(iC);
        ALL.cellDepth(cc) = S.usr.tetrodeDepths(iC);
        ALL.ID(cc) = s_out.cq.id(iC); ALL.Lr(cc) = s_out.cq.lr(iC); ALL.ampl(cc) = s_out.cq.ampl(iC);
         
        ALL.bigP(cc,:) = P; F_acf = F; % spike spectrum
        ALL.sessno(cc) = iS;
        ALL.ratID(cc) = fd_extra.ratID_num(iS);
        
        % shuffles for spike spectrum
        clear this_shufP;
        for iShuf = cfg_master.nShuf:-1:1
            
            % shuffle spike train
            spk_shuf = rand(length(spk_t), 1);
            spk_shuf = tbin_edges(1) + spk_shuf .* (tbin_edges(end) - tbin_edges(1));
            
            % re-trialify
            clear data;
            for iT = 1:length(trial_starts) % this bit should be made into a function
                this_spk = spk_shuf - trial_starts(iT);
                keep = this_spk >= 0 & this_spk < cfg.trial_len;
                data(iT).times = this_spk(keep);
            end
            
            [P,F,R] = mtspectrumpt(data, cfg.params);
            P = nanmean(P, 2) ./ length(spk_shuf);
            keep = F > 1;
            P = 10*log10(P(keep)); F = F(keep);
            
            this_shufP(iShuf,:) = P;
        end % of shuffles
        
        ALL.bigP_shufmean(cc,:) = nanmean(this_shufP);
        ALL.bigP_shufSD(cc,:) = nanstd(this_shufP);
        
        if cfg_master.plot
            hold on;
            plot(F, ALL.bigP_shufmean(cc,:), 'r', 'LineWidth', 2);
            plot(F, ALL.bigP_shufmean(cc,:) + ALL.bigP_shufSD(cc,:), 'r');
            plot(F, ALL.bigP_shufmean(cc,:) - ALL.bigP_shufSD(cc,:), 'r');
            drawnow;
        end
        
        % get STA by binning on CSC timebase
        this_spk_binned = zeros(size(csc.tvec));
        idx = nearest_idx3(spk_t, csc.tvec);
        this_spk_binned(idx) = 1;
        [xc, tvec] = xcorr(csc.data, this_spk_binned, 1000);
        tvec = tvec .* median(diff(csc.tvec));
        
        ALL.STA(cc,:) = xc;
        
        if cfg_master.plot
            subplot(222);
            plot(tvec, xc);
            title('raw STA');
        end
        
        % STA spectrum
        this_Fs = 1 ./ median(diff(csc.tvec));
        [P,F] = pwelch(xc,length(xc), [], [], this_Fs);
        
        if cfg_master.plot
            subplot(224);
            plot(F, 10*log10(P)); xlim([0 100]);
            title('STA spectrum');
        end
        
        ALL.STAp(cc,:) = P;
        
        % do STA shuffles
        clear this_shufSTA this_shufSTAp;
        for iShuf = cfg.nShuf:-1:1
            spk_binned_shuf = this_spk_binned(randperm(length(this_spk_binned)));
            this_shufSTA(iShuf,:) = xcorr(csc.data, spk_binned_shuf,1000);
            this_shufSTAp(iShuf,:) = pwelch(this_shufSTA(iShuf,:), length(this_shufSTA(iShuf,:)), [], [], this_Fs);
        end
        ALL.STA_shufmean(cc,:) = nanmean(this_shufSTA); ALL.STAp_shufmean(cc,:) = nanmean(this_shufSTAp);
        ALL.STA_shufSD(cc,:) = nanstd(this_shufSTA); ALL.STAp_shufSD(cc,:) = nanstd(this_shufSTAp);
        
        if cfg_master.plot
            subplot(222); hold on;
            plot(tvec,ALL.STA_shufmean(cc,:),'r');
            plot(tvec,ALL.STA_shufmean(cc,:) + 2*ALL.STA_shufSD(cc,:), 'r:');
            plot(tvec,ALL.STA_shufmean(cc,:) - 2*ALL.STA_shufSD(cc,:), 'r:');
            
            subplot(224); hold on;
            plot(F,10*log10(ALL.STAp_shufmean(cc,:)), 'r');
            plot(F,10*log10(ALL.STAp_shufmean(cc,:) + 2*ALL.STAp_shufSD(cc,:)), 'r:');
            plot(F,10*log10(ALL.STAp_shufmean(cc,:) - 2*ALL.STAp_shufSD(cc,:)), 'r:');
            
            drawnow;
        end
        
        cc = cc + 1; % advance cell count
    end % of cell loop
    
end % of session loop

%% compute some z-scores and plot
f_int = 1:0.5:100;

bigPi = interp1(F_acf, ALL.bigP', f_int);
bigPz = (ALL.bigP - ALL.bigP_shufmean) ./ ALL.bigP_shufSD;
bigPzi = interp1(F_acf, bigPz', f_int);

STApi = interp1(F, ALL.STAp', f_int);
STApz = (ALL.STAp - ALL.STAp_shufmean) ./ ALL.STAp_shufSD;
STApzi = interp1(F, STApz', f_int);

keep = ~isnan(sum(bigPz, 2));
STApi = STApi(:, keep); bigPi = bigPi(:, keep);
STApzi = STApzi(:, keep); bigPzi = bigPzi(:, keep);

%%
figure;
subplot(221);
imagesc(f_int, 1:size(bigPz,1), bigPzi'); caxis([-10 10]);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100); grid on;
title('spike spectrum (z)'); xlabel('frequency (Hz)'); ylabel('cell #');

subplot(222);
imagesc(f_int, 1:size(STApz,1), STApzi'); caxis([-10 10]);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100); grid on;
title('STA spectrum (z)'); xlabel('frequency (Hz)'); ylabel('cell #');

subplot(223);
ss_count = sum(bigPzi > 1.96, 2);
plot(f_int, ss_count ./ size(bigPzi, 2), 'k', 'LineWidth', 2);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100); grid on; box off;
xlabel('frequency (Hz)'); ylabel('%cells significant');

subplot(224);
sta_count = sum(STApzi > 1.96, 2);
plot(f_int, sta_count ./ size(bigPzi, 2), 'k', 'LineWidth', 2);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100); grid on; box off;
xlabel('frequency (Hz)'); ylabel('%cells significant');

%%
figure;

subplot(221);
cc = corrcoef(bigPi');
imagesc(f_int, f_int, cc);
caxis([0 0.5]);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100, 'YLim', [1 100], 'YTick', 0:10:100); grid on;
title('spike spectrum (raw) xcorr'); colorbar;

subplot(222);
cc = corrcoef(STApi');
imagesc(f_int, f_int, cc);
caxis([-0.25 1]);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100, 'YLim', [1 100], 'YTick', 0:10:100); grid on;
title('STA spectrum (raw) xcorr'); colorbar;

subplot(223);
cc = corrcoef(bigPzi');
imagesc(f_int, f_int, cc);
caxis([0 0.5]);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100, 'YLim', [1 100], 'YTick', 0:10:100); grid on;
title('spike spectrum (z) xcorr'); colorbar;

subplot(224);
cc = corrcoef(STApzi');
imagesc(f_int, f_int, cc);
caxis([-0.25 1]);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100, 'YLim', [1 100], 'YTick', 0:10:100); grid on;
title('STA spectrum (z) xcorr'); colorbar;

%%
clear cc cc_raw;
for iF = 1:length(f_int)
    
    [temp1,temp2] = corrcoef(STApi(iF,:),bigPi(iF,:));
    cc_raw(iF) = temp1(1,2);
    pval_raw(iF) = temp2(1,2);
    
    [temp1,temp2] = corrcoef(STApzi(iF,:),bigPzi(iF,:));
    cc(iF) = temp1(1,2);
    pval(iF) = temp2(1,2);
    
end

figure;

subplot(211);
plot(f_int, cc, 'k', 'LineWidth', 2);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100); grid on;
xlabel('frequency (Hz)'); ylabel('correlation'); title('SS x STA (z)');

xval = f_int;
keep = pval < 0.05;
hold on; plot(xval(keep),cc(keep),'.b','MarkerSize',20);

subplot(212);
plot(f_int, cc_raw, 'k', 'LineWidth', 2);
set(gca, 'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out', 'XLim', [1 100], 'XTick', 0:10:100); grid on;
xlabel('frequency (Hz)'); ylabel('correlation'); title('SS x STA (raw)');

xval = f_int;
keep = pval_raw < 0.05;
hold on; plot(xval(keep),cc_raw(keep),'.b','MarkerSize',20);

%%
function [acf,tvec] = ComputeACF(cfg,spk_binned)
if isfield(cfg,'maxlag')
    [acf,tvec] = xcorr(spk_binned,spk_binned,cfg.maxlag);
else
    [acf,tvec] = xcorr(spk_binned,spk_binned);
end
tvec = tvec.*cfg.binsize;
acf(ceil(length(acf)/2)) = 0;
end