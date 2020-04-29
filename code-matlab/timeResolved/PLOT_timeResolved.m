%% quick plot for single cells
tfr_baseline = [-5 -4.5];
for iC = 1:cc

    figure;
    
    % session-specgram
    if isfield(ALL, 'sessionTFR')
        subplot(331);
        t_idx = nearest_idx3(tfr_baseline, ALL.sessionTFR(iC).time);
        baseline = sq(nanmean(ALL.sessionTFR(iC).powspctrm(:, :, t_idx(1):t_idx(2)), 3));
        baseline = repmat(baseline', [1 size(ALL.sessionTFR(iC).powspctrm, 3)]);
        
        imagesc(ALL.sessionTFR(iC).time, ALL.sessionTFR(iC).freq, sq(ALL.sessionTFR(iC).powspctrm)-baseline); axis xy; colorbar;
        set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('time'); ylabel('frequency (Hz)'); title('LFP spectrogram');
    end
    
    % spike spectrum -- COULD SMOOTH THIS SOME
    subplot(332);
    plot(ALL.spkSpec_freq, medfilt1(ALL.spkSpec(iC, :),11,'truncate'), 'k', 'LineWidth', 2);
    hold on;
    plot(ALL.spkSpec_freq, medfilt1(ALL.spkSpec_shufmean(iC, :),11) + medfilt1(ALL.spkSpec_shufSD(iC, :),11), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    plot(ALL.spkSpec_freq, medfilt1(ALL.spkSpec_shufmean(iC, :),11) - medfilt1(ALL.spkSpec_shufSD(iC, :),11), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('frequency (Hz)'); ylabel('power'); title('spike spectrum');
    
    % STA (raw)
    subplot(333);
    plot(ALL.STAtime, ALL.STA(iC, :), 'k', 'LineWidth', 1);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-0.5 0.5],'XTick',-0.5:0.1:0.5); xlabel('time (s)'); ylabel('LFP'); title('STA');
    set(gca,'XTickLabel',{'-0.5','','','','','0','','','','','0.5'});
    grid on;
    
    subplot(334);
    plot(ALL.ppc_freq, ALL.ppc(iC, :), 'k', 'LineWidth', 2);
    hold on;
    plot(ALL.ppc_freq, ALL.ppc_shufmean(iC, :) + ALL.ppc_shufSD(iC, :), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    plot(ALL.ppc_freq, ALL.ppc_shufmean(iC, :) - ALL.ppc_shufSD(iC, :), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('frequency (Hz)'); ylabel('ppc'); title('spike-field ppc');
    % try to recover this average from ppc-gram and spike count?
    
    subplot(335);
    %imagesc(ALL.trPPCtime, ALL.trPPCfreq, sq(ALL.trPPC(iC,:,:))); axis xy; colorbar;
    imagesc(ALL.trPPCtime, ALL.trPPCfreq, (sq(ALL.trPPC(iC,:,:)) - sq(ALL.trPPC_shufmean(iC,:,:))) ./ sq(ALL.trPPC_shufSD(iC,:,:))); axis xy; colorbar;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('frequency (Hz)'); title('PPC (z)');
    %caxis([0 8*max(ALL.ppc(iC,:))])
    caxis(caxis ./ 2);
    
    subplot(336);
    plot(ALL.STAp_freq,ALL.STAp(iC,:), 'k', 'LineWidth', 2);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('frequency (Hz)'); ylabel('power'); title('STA power');
    
    subplot(337); bar(ALL.trPPCtime, sq(ALL.trPPCn(iC,1,:))); box off;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('count'); title('spike count histogram');
    
    subplot(338)
    imagesc(ALL.tr_time-5.5, ALL.tr_freq, sq(ALL.trP(iC,:,:))'); axis xy; colorbar;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('frequency (Hz)'); title('spike spectrum (raw)');
    % how to normalize this by spike count? need histogram on chronux timebase...
    
    this_trp = sq(ALL.trP(iC,:,:))';
    sc = repmat(ALL.trP_hist(iC,:), [size(this_trp,1) 1]);
    this_trp = this_trp ./ sc;
    
    subplot(339)
    imagesc(ALL.tr_time-5.5, ALL.tr_freq, this_trp); axis xy; colorbar;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('frequency (Hz)'); title('spike spectrum (n)');
    
    drawnow;
    %pause; clf;

end


%% some averages
what = {'FSI', 'MSNonly'};

ib = 1:0.5:100; % new frequency basis
hh = 0.1; % histogram height
ph = 1; % number of vertical subplots (controls height)
c_scale = 0.25;

for iW = 1:length(what)
    
    switch what{iW} % select appropriate cell idxs
        case 'FSI'
            keep_idx = find(ALL.cellType == 2);
        case 'MSNonly'
            keep_idx = find(ALL.cellType == 1);
        case 'nonFSI'
            keep_idx = find(ALL.cellType ~= 2);
        case 'all'
            keep_idx = 1:length(ALL.cellType);
    end
    nCells = length(keep_idx);
    fprintf('%d cells of %d (%.2f %%) are %s\n', nCells, length(ALL.cellType), (nCells ./ length(ALL.cellType)) .* 100, what{iW});
    
    figure;
    
    % spike spectra
    if isfield(ALL, 'spkSpec')
        this_ss = ALL.spkSpec(keep_idx, :)';
        this_ss = interp1(ALL.spkSpec_freq, this_ss, ib, 'linear')';
        
        this_ssN = (ALL.spkSpec(keep_idx, :) - ALL.spkSpec_shufmean(keep_idx, :)) ./ ALL.spkSpec_shufSD(keep_idx, :);
        this_ssN = interp1(ALL.spkSpec_freq, this_ssN', ib, 'linear')';
        
        as = subplot(ph, 3, 1);
        ax1 = axes('Position', [as.Position(1) as.Position(2) + hh * as.Position(4) as.Position(3) (1-hh) * as.Position(4)]);
        ax2 = axes('Position', [as.Position(1) as.Position(2) as.Position(3) hh * as.Position(4)]);
        delete(as);
        
        axes(ax1); imagesc(ib, 1:nCells, this_ss); axis xy
        set(gca, 'LineWidth', 1, 'XLim', [1 100], 'XTick', [1 10:10:100], 'XTickLabel', [], 'FontSize', 18, 'YTick', []);
        ylabel('cell #'); title(sprintf('spike spectrum (raw, %s)', what{iW})); grid on;
        %caxis(caxis.*c_scale);
        caxis([0 1*10e-3]);
        
        axes(ax2); plot(ib, nanmean(this_ss), 'k'); box off;
        set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
            'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
            'FontSize', 18, 'YTick', []);
        grid on; xlabel('frequency (Hz)');
        ylim([0.25e-3 2.5e-3]);
    end
    
    % ppc
    this_ppc = ALL.ppc(keep_idx, :)';
    this_ppc = interp1(ALL.ppc_freq, this_ppc, ib, 'linear')';
    
    as = subplot(ph, 3, 2);
    ax1 = axes('Position', [as.Position(1) as.Position(2) + hh * as.Position(4) as.Position(3) (1-hh) * as.Position(4)]); 
    ax2 = axes('Position', [as.Position(1) as.Position(2) as.Position(3) hh * as.Position(4)]);
    delete(as);
    
    axes(ax1); imagesc(ib, 1:nCells, this_ppc); axis xy
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'XTick', [1 10:10:100], 'XTickLabel', [], 'FontSize', 18, 'YTick', []);
    ylabel('cell #'); title(sprintf('ppc (raw, %s)', what{iW})); grid on;
    %caxis(caxis.*c_scale);
    caxis([0 0.04]);
    
    axes(ax2); plot(ib, nanmean(this_ppc), 'k'); box off;
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', []);
    grid on; xlabel('frequency (Hz)');
    ylim([-0.001 0.01]);
    
    % stap
    if isfield(ALL, 'STAp')
    this_sta = ALL.STAp(keep_idx, :)';
    this_sta = interp1(ALL.STAp_freq, this_sta, ib, 'linear')';
    
    as = subplot(ph, 3, 3);
    ax1 = axes('Position', [as.Position(1) as.Position(2) + hh * as.Position(4) as.Position(3) (1-hh) * as.Position(4)]); 
    ax2 = axes('Position', [as.Position(1) as.Position(2) as.Position(3) hh * as.Position(4)]);
    delete(as);
    
    axes(ax1); imagesc(ib, 1:nCells, this_sta); axis xy
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'XTick', [1 10:10:100], 'XTickLabel', [], 'FontSize', 18, 'YTick', []);
    ylabel('cell #'); title(sprintf('STA power (raw, %s)', what{iW})); grid on;
    %caxis(caxis.*c_scale);
    caxis([0 10]);
    
    axes(ax2); plot(ib, nanmean(this_sta), 'k'); box off;
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', []);
    grid on; xlabel('frequency (Hz)'); 
    ylim([0 2]);
    end
    
    figure; % z-scored
    
    if isfield(ALL, 'spkSpec')
    as = subplot(ph, 3, 1);
    ax1 = axes('Position', [as.Position(1) as.Position(2) + hh * as.Position(4) as.Position(3) (1-hh) * as.Position(4)]); 
    ax2 = axes('Position', [as.Position(1) as.Position(2) as.Position(3) hh * as.Position(4)]);
    delete(as);
    
    axes(ax1); imagesc(ib, 1:nCells, this_ssN); axis xy
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'XTick', [1 10:10:100], 'XTickLabel', [], 'FontSize', 18, 'YTick', []);
    ylabel('cell #'); title(sprintf('spike spectrum (z, %s)', what{iW})); grid on;
    caxis(caxis.*c_scale);
    
    axes(ax2); plot(ib, nanmean(this_ssN), 'k'); box off;
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', []);
    grid on; xlabel('frequency (Hz)'); 
    end
    
    if isfield(ALL, 'ppc_shufmean')
    this_ppcN = (ALL.ppc(keep_idx, :) - ALL.ppc_shufmean(keep_idx, :)) ./ ALL.ppc_shufSD(keep_idx, :);
    this_ppcN = interp1(ALL.ppc_freq, this_ppcN', ib, 'linear')';
        
    as = subplot(ph, 3, 2);
    ax1 = axes('Position', [as.Position(1) as.Position(2) + 0.1 * as.Position(4) as.Position(3) 0.9 * as.Position(4)]); 
    ax2 = axes('Position', [as.Position(1) as.Position(2) as.Position(3) 0.1 * as.Position(4)]);
    delete(as);
    
    axes(ax1); imagesc(ib, 1:nCells, this_ppcN); axis xy
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'XTick', [1 10:10:100], 'XTickLabel', [], 'FontSize', 18, 'YTick', []);
    ylabel('cell #'); title(sprintf('ppc (z, %s)', what{iW})); grid on;
    caxis(caxis.*c_scale);
    
    axes(ax2); plot(ib, nanmean(this_ppcN), 'k'); box off;
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', []);
    grid on; xlabel('frequency (Hz)'); 
    end
    
    figure; % histograms/counts
    
    
    if isfield(ALL, 'spkSpec')
        subplot(331)
        hb = bar(ib, nanmean(this_ssN > 1.96), 'k'); box off; set(hb, 'EdgeColor', 'none');
        set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
            'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
            'FontSize', 18, 'YLim', [0 1], 'YTick', 0:0.25:1, 'YTickLabel', {'0','','0.5','','1'});
        xlabel('frequency (Hz)'); grid on; title('SS prop. significant');
    end
    
    if isfield(ALL, 'ppc_shufmean')
        subplot(332)
        hb = bar(ib, nanmean(this_ppcN > 1.96), 'k'); box off; set(hb, 'EdgeColor', 'none');
        set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
            'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
            'FontSize', 18, 'YLim', [0 1], 'YTick', 0:0.25:1, 'YTickLabel', {'0','','0.5','','1'});
        
        xlabel('frequency (Hz)'); grid on; title('PPC prop. significant');
    end
      
    % phases -- need to get histogram for each frequency (PPC version)
    this_phi = ALL.ppc_ang(keep_idx, :)';
    this_phi = interp1(ALL.ppc_freq, this_phi, ib, 'nearest')';
    
    this_STAa = ALL.STAa(keep_idx, :)';
    this_STAa = interp1(ALL.STAp_freq, this_STAa, ib, 'linear')';
    
    this_r = ALL.ppc_r(keep_idx, :)';
    this_r = interp1(ALL.ppc_freq, this_r, ib, 'linear')';
    
    phi_bin_edges = -pi:pi/3:pi; phi_bin_centers = phi_bin_edges(1:end-1) + mode(diff(phi_bin_edges))/2;
    clear phi_hist phi_hist_STA;
    for iB = size(this_phi, 2):-1:1
        phi_hist(iB, :) = histc(this_phi(this_ppcN(:, iB) > 1.96, iB), phi_bin_edges);
        phi_hist_STA(iB, :) = histc(this_STAa(this_ppcN(:, iB) > 1.96, iB), phi_bin_edges);
    end
    phi_hist = phi_hist(:, 1:end-1); phi_hist_STA = phi_hist_STA(:, 1:end-1);
    
    subplot(334)
    imagesc(ib, phi_bin_centers, phi_hist');
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YLim', [-pi pi], 'YTick', [phi_bin_centers(1) phi_bin_centers(end)], 'YTickLabel', {'-pi', 'pi'});
    xlabel('frequency (Hz)'); grid on; title('raw phase histogram (PPC)'); axis tight;
    
    phi_histN = phi_hist ./ repmat(nansum(phi_hist, 2), [1 size(phi_hist, 2)]);
    subplot(335)
    imagesc(ib,pb,phi_histN');
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YLim', [-pi pi], 'YTick', [phi_bin_centers(1) phi_bin_centers(end)], 'YTickLabel', {'-pi', 'pi'});
    xlabel('frequency (Hz)'); grid on; title('normalized phase histogram (PPC)'); axis tight;
    
    subplot(337)
    imagesc(ib,pb,phi_hist_STA');
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YLim', [-pi pi], 'YTick', [phi_bin_centers(1) phi_bin_centers(end)], 'YTickLabel', {'-pi', 'pi'});
    xlabel('frequency (Hz)'); grid on; title('raw phase histogram (STA)'); axis tight;
    
    phi_histN = phi_hist_STA ./ repmat(nansum(phi_hist_STA, 2), [1 size(phi_hist_STA, 2)]);
    subplot(338)
    imagesc(ib,pb,phi_histN');
    set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YLim', [-pi pi], 'YTick', [phi_bin_centers(1) phi_bin_centers(end)], 'YTickLabel', {'-pi', 'pi'});
    xlabel('frequency (Hz)'); grid on; title('normalized phase histogram (STA)'); axis tight;

    % phases -- need to get histogram for each frequency
%     if isfield(ALL, 'ppc_ang')
%         
%         this_phi = ALL.ppc_ang(keep_idx, :)';
%         this_phi = interp1(ALL.ppc_freq, this_phi, ib, 'nearest')';
%         for iB = length(ib):-1:1
%             [phi_hist(iB, :), pb] = hist(this_phi(this_ppcN(:, iB) > 1.96, iB), 6);
%         end
%         
%         subplot(334)
%         imagesc(ib,pb,phi_hist');
%         set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
%             'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
%             'FontSize', 18, 'YLim', [-pi pi], 'YTick', [-pi pi], 'YTickLabel', {'-pi', 'pi'});
%         xlabel('frequency (Hz)'); grid on; title('raw phase histogram');
%         
%         phi_histN = phi_hist ./ repmat(nansum(phi_hist, 2), [1 size(phi_hist, 2)]);
%         subplot(335)
%         imagesc(ib,pb,phi_histN');
%         set(gca, 'LineWidth', 1, 'XLim', [1 100], 'TickDir', 'out', 'XTick', [1 10:10:100], ...
%             'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
%             'FontSize', 18, 'YLim', [-pi pi], 'YTick', [-pi pi], 'YTickLabel', {'-pi', 'pi'});
%         xlabel('frequency (Hz)'); grid on; title('normalized phase histogram');
%         
%     end % of phase histogram

    
    figure; % freq specific histograms (raw)
    
    fb = {[2.5 5.5], [6.5 9.5], [13.5 25.5], [39.5 65.5], [65.5 90.5]};
    fn = {'delta', 'theta', 'beta', 'lgamma', 'hgamma'};
    nBins = 100;
    
    for iF = 1:length(fb)
             
        this_ib = find(ib >= fb{iF}(1) & ib < fb{iF}(2));
        
        % ss
        if isfield(ALL, 'spkSpec')
            
            temp = nanmean(this_ssN(:, this_ib), 2);
            sig_idx = temp > 1.96;
            
            if sum(sig_idx) > 0
                
                ss_avg = nanmean(this_ss(sig_idx, this_ib), 2); % could consider doing max...
                bin_edgesLog = linspace(floor(10*log10(min(ss_avg)))/10, ceil(10*log10(max(ss_avg)))/10, nBins);
                bin_centersLog = bin_edgesLog(1:end-1) + mode(diff(bin_edgesLog)) ./ 2;
                this_h = histc(log10(ss_avg), bin_edgesLog); this_h = this_h(1:end-1);
                
                subplot(4, 5, iF)
                hb = bar(bin_centersLog, this_h, 'k'); box off; set(hb, 'EdgeColor', 'none');
                set(gca, 'LineWidth', 1, 'TickDir', 'out', ...
                    'FontSize', 18), ...
                    xlabel('log power'); title(sprintf('%s ss (n = %d, %.2f%%)', fn{iF}, sum(sig_idx), (sum(sig_idx) ./ nCells) .* 100));
                
                ss_avg = nanmean(this_ssN(sig_idx, this_ib), 2); % could consider doing max...
                %[this_h, this_b] = hist(ss_avg, nBins);
                bin_edgesLog = linspace(floor(10*log10(min(ss_avg)))/10, ceil(10*log10(max(ss_avg)))/10, nBins);
                bin_centersLog = bin_edgesLog(1:end-1) + mode(diff(bin_edgesLog)) ./ 2;
                this_h = histc(log10(ss_avg), bin_edgesLog); this_h = this_h(1:end-1);
                
                subplot(4, 5, iF+5)
                hb = bar(bin_centersLog, this_h, 'k'); box off; set(hb, 'EdgeColor', 'none');
                set(gca, 'LineWidth', 1, 'TickDir', 'out', ...
                    'FontSize', 18), ...
                    xlabel('log power (z)'); title(sprintf('%s ss (n = %d, %.2f%%)', fn{iF}, sum(sig_idx), (sum(sig_idx) ./ nCells) .* 100));
                
            end
        end
     
        % ppc
    
        temp = nanmean(this_ppcN(:, this_ib), 2);
        sig_idx = temp > 1.96;
    
        if sum(sig_idx) > 0
        
            ss_avg = nanmean(this_ppc(sig_idx, this_ib), 2); % could consider doing max...
            bin_edgesLog = linspace(floor(10*log10(min(ss_avg)))/10, ceil(10*log10(max(ss_avg)))/10, nBins);
            bin_centersLog = bin_edgesLog(1:end-1) + mode(diff(bin_edgesLog)) ./ 2;
            this_h = histc(log10(ss_avg), bin_edgesLog); this_h = this_h(1:end-1);
            
            subplot(4, 5, iF+10)
            hb = bar(bin_centersLog, this_h, 'k'); box off; set(hb, 'EdgeColor', 'none');
            set(gca, 'LineWidth', 1, 'TickDir', 'out', ...
                'FontSize', 18), ...
                xlabel('log PPC'); title(sprintf('%s ppc (n = %d, %.2f%%)', fn{iF}, sum(sig_idx), (sum(sig_idx) ./ nCells) .* 100));
                          
            ss_avg = nanmean(this_ppcN(sig_idx, this_ib), 2); % could consider doing max...
            bin_edgesLog = linspace(floor(10*log10(min(ss_avg)))/10, ceil(10*log10(max(ss_avg)))/10, nBins);
            bin_centersLog = bin_edgesLog(1:end-1) + mode(diff(bin_edgesLog)) ./ 2;
            this_h = histc(log10(ss_avg), bin_edgesLog); this_h = this_h(1:end-1);
            
            subplot(4, 5, iF+15)
            hb = bar(bin_centersLog, this_h, 'k'); box off; set(hb, 'EdgeColor', 'none');
            set(gca, 'LineWidth', 1, 'TickDir', 'out', ...
                'FontSize', 18), ...
                xlabel('log PPC (z)'); title(sprintf('%s ppcN (n = %d, %.2f%%)', fn{iF}, sum(sig_idx), (sum(sig_idx) ./ nCells) .* 100));
               
        end
    
    end
     
    % phase histograms using PPC measures

    if isfield(ALL, 'ppc_ang')
    figure;
    
    for iF = 1:length(fb)
        
        this_ib = find(ib >= fb{iF}(1) & ib < fb{iF}(2));
        
        % find significant PPC cells
        temp = nanmean(this_ppcN(:, this_ib), 2);
        %temp = max(this_ppcN(:, this_ib), [], 2);
        sig_idx = temp > 1.96; % idx of significant cells
        
        % find angle for largest ppc within band
        %this_ppc_band = nanmean(this_ppc(:, this_ib), 2);
        [this_ppc_band, temp_idx] = max(this_ppc(sig_idx, this_ib)');
        this_phi_mat = this_phi(sig_idx, :);
        clear this_phi_band
        for iC = size(this_phi_mat, 1):-1:1 % get phi at freq with largest PPC
            this_phi_band(iC) = this_phi_mat(iC, temp_idx(iC));
        end
        
        % polar plot with length equal to ppc
        subplot(4, 5, iF);
        for iC = 1:length(this_ppc_band)
            %h = polarplot([0 this_phi_band(iC)], [0 log(this_ppc_band(iC))], '-k'); set(h, 'LineWidth', 1);
            h = polarplot([0 this_phi_band(iC)], [0 this_ppc_band(iC)], '-k'); set(h, 'LineWidth', 1);
            hold on;
        end
        %set(gca, 'FontSize', 18, 'RLim', [0 0.1], 'RTickLabel', {}, 'ThetaTickLabel', {}); 
        set(gca, 'FontSize', 18, 'RTickLabel', {}, 'ThetaTickLabel', {});
        title(sprintf('%s (n = %d) raw PPC', fn{iF}, length(this_phi_band)));
        
        % regular histograms (not scaled by PPC)
        subplot(4, 5, iF+5);
        phi_bin_edges = -pi:pi/6:pi; phi_bin_centers = phi_bin_edges(1:end-1) + mode(diff(phi_bin_edges))/2;
        this_hist = histc(this_phi_band, phi_bin_edges); this_hist = this_hist(1:end-1);
        hb = bar(phi_bin_centers, this_hist, 'k'); box off; set(hb, 'EdgeColor', 'none'); set(hb,'BarWidth', 1)
        set(gca, 'XTick', phi_bin_edges, 'XTickLabel', []); axis tight;
        set(gca, 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 18);
        xlabel('phi');
                       
        
        % z-scored PPC
        subplot(4, 5, iF+10);
        [this_ppc_band, temp_idx] = max(this_ppcN(sig_idx, this_ib)');
        this_phi_mat = this_phi(sig_idx, :);
        clear this_phi_band
        for iC = size(this_phi_mat, 1):-1:1 % get phi at freq with largest PPC
            this_phi_band(iC) = this_phi_mat(iC, temp_idx(iC));
        end

        for iC = 1:length(this_ppc_band)
            %h = polarplot([0 this_phi_band(iC)], [0 log(this_ppc_band(iC))], '-k'); set(h, 'LineWidth', 1);
            h = polarplot([0 this_phi_band(iC)], [0 this_ppc_band(iC)], '-k'); set(h, 'LineWidth', 1);
            hold on;
        end
        %set(gca, 'FontSize', 18, 'RLim', [0 0.1], 'RTickLabel', {}, 'ThetaTickLabel', {}); 
        set(gca, 'FontSize', 18, 'RTickLabel', {}, 'ThetaTickLabel', {});
        title(sprintf('%s (n = %d) PPCz', fn{iF}, length(this_phi_band)));
        
        % regular histograms (not shown by PPC)
        subplot(4, 5, iF+15);
        phi_bin_edges = -pi:pi/6:pi; phi_bin_centers = phi_bin_edges(1:end-1) + mode(diff(phi_bin_edges))/2;
        this_hist = histc(this_phi_band, phi_bin_edges); this_hist = this_hist(1:end-1);
        hb = bar(phi_bin_centers, this_hist, 'k'); box off; set(hb, 'EdgeColor', 'none'); set(hb,'BarWidth', 1)
        set(gca, 'XTick', phi_bin_edges, 'XTickLabel', []); axis tight;
        set(gca, 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 18);
        xlabel('phi');
        
    end % of freq bands
    end
    
    % phase histograms using STA and r   
    figure;
    
    for iF = 1:length(fb)
        
        this_ib = find(ib >= fb{iF}(1) & ib < fb{iF}(2));
        
        % find significant PPC cells
        temp = nanmean(this_ppcN(:, this_ib), 2);
        %temp = max(this_ppcN(:, this_ib), [], 2);
        sig_idx = temp > 1.96; % idx of significant cells
        
        % find angle for largest ppc within band
        %this_ppc_band = nanmean(this_ppc(:, this_ib), 2);
        [this_ppc_band, temp_idx] = max(this_ppc(sig_idx, this_ib)');
        this_phi_mat = this_STAa(sig_idx, :);
        clear this_phi_band
        for iC = size(this_phi_mat, 1):-1:1 % get phi at freq with largest PPC
            this_phi_band(iC) = this_phi_mat(iC, temp_idx(iC));
        end
        
        % polar plot with length equal to ppc
        subplot(4, 5, iF);
        for iC = 1:length(this_ppc_band)
            %h = polarplot([0 this_phi_band(iC)], [0 log(this_ppc_band(iC))], '-k'); set(h, 'LineWidth', 1);
            h = polarplot([0 this_phi_band(iC)], [0 this_ppc_band(iC)], '-k'); set(h, 'LineWidth', 1);
            hold on;
        end
        %set(gca, 'FontSize', 18, 'RLim', [0 0.1], 'RTickLabel', {}, 'ThetaTickLabel', {});
        set(gca, 'FontSize', 18, 'RTickLabel', {}, 'ThetaTickLabel', {});
        title(sprintf('%s (n = %d) raw |r|', fn{iF}, length(this_phi_band)));
        
        % regular histograms (not shown by PPC)
        subplot(4, 5, iF+5);
        phi_bin_edges = -pi:pi/6:pi; phi_bin_centers = phi_bin_edges(1:end-1) + mode(diff(phi_bin_edges))/2;
        this_hist = histc(this_phi_band, phi_bin_edges); this_hist = this_hist(1:end-1);
        hb = bar(phi_bin_centers, this_hist, 'k'); box off; set(hb, 'EdgeColor', 'none'); set(hb,'BarWidth', 1)
        set(gca, 'XTick', phi_bin_edges, 'XTickLabel', []); axis tight;
        set(gca, 'LineWidth', 1, 'TickDir', 'out', 'FontSize', 18);
        xlabel('phi');
        
    end % of freq bands
    
    
end % of cell types


%% correlations
ca = [-0.1 0.9];

figure;

% ss v ss
subplot(331)
cc = corrcoef(this_ss);
imagesc(ib, ib, cc); axis xy;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', [1 10:10:100], ...
        'YTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'});
grid on; ylabel('SS (raw)'); title(sprintf('SS (raw)')); caxis(ca);
    
% ppc v ss
subplot(332)
cc = corrcoef([this_ss this_ppc]);
%cc = cc(1:length(ib), length(ib)+1:end);
cc = cc(length(ib)+1:end, 1:length(ib));
imagesc(ib, ib, cc); axis xy;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', [1 10:10:100], ...
        'YTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'});
    grid on; title(sprintf('PPC (raw)')); caxis(ca);
    
subplot(334) % just the diagonal
plot(cc(logical(eye(size(cc, 1)))), 'k'); box off;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], 'XLim', [1 100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YLim', [0 1]); ylabel('PPC correlation'); grid on;
    
% sta v ss
subplot(333)
cc = corrcoef([this_ss this_sta]);
%cc = cc(1:length(ib), length(ib)+1:end);
cc = cc(length(ib)+1:end, 1:length(ib));
imagesc(ib, ib, cc); axis xy;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', [1 10:10:100], ...
        'YTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'});
    grid on; title(sprintf('STA (raw)')); caxis(ca);
     
subplot(337) % just the diagonal
plot(cc(logical(eye(size(cc, 1)))), 'k'); box off;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], 'XLim', [1 100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YLim', [0 1]); ylabel('STA correlation'); grid on;
    
% ppc v ppc
subplot(335)
cc = corrcoef(this_ppc);
imagesc(ib, ib, cc); axis xy;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', [1 10:10:100], ...
        'YTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'});
grid on; caxis(ca);

% ppc v sta
subplot(336)
cc = corrcoef([this_ppc this_sta]);
%cc = cc(1:length(ib), length(ib)+1:end);
cc = cc(length(ib)+1:end, 1:length(ib));
imagesc(ib, ib, cc); axis xy;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', [1 10:10:100], ...
        'YTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'});
grid on; caxis(ca);
    
subplot(338) % just the diagonal
plot(cc(logical(eye(size(cc, 1)))), 'k'); box off;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], 'XLim', [1 100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YLim', [0 1]); grid on;

% sta v sta
subplot(339)
cc = corrcoef(this_sta);
imagesc(ib, ib, cc); axis xy;
set(gca, 'LineWidth', 1, 'TickDir', 'out', 'XTick', [1 10:10:100], ...
        'XTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'}, ...
        'FontSize', 18, 'YTick', [1 10:10:100], ...
        'YTickLabel', {'1', '', '20', '', '40', '', '60', '', '80', '', '100'});
grid on; caxis(ca);

