cd('D:\RandomVstrAnalysis\ft_results');
% cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
hg = [5,100];
pk_thresh = -1;
num_control_splits = 100;

load('./significant_dif.mat');
%%
%Generate FSI_stuff
for iC = 1:length(sig_fsi_labels)
   %generate correct file to read
   this_label = sig_fsi_labels{iC};
   fname = split(sig_fsi_labels{iC},'-');
   fname = join(fname(1:4),'-');
   fname = cat(2, fname{1}, '_ft_spec.mat');
   load(fname);
   fsi_labels  = od.label(od.cell_type == 2);
   for iF = 1:length(fsi_labels)
      if matches(fsi_labels(iF), this_label)
        close all;
        fig = figure('WindowState', 'maximized');
        % Plot STA
        ax1 = subplot(3,1,1);
        hold on;
        p1 = plot(ax1, od.fsi_res.near_spec{iF}.sta_time, ...
            od.fsi_res.near_lfr_spec{iF}.sta_vals, 'Color', 'red');
        p2 = plot(ax1, od.fsi_res.near_spec{iF}.sta_time, ...
            od.fsi_res.near_hfr_spec{iF}.sta_vals, 'Color', 'green');
        ax1.Box = 'off';
        ax1.YTick = [];
        ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
        ax1.Title.String = 'Spike Triggered Average (STA)';
        ax1.XLabel.String = 'Time (in sec)';
        ax1.Title.FontSize = 20;
        ax1.XAxis.FontSize = 18;
        leg1 = legend({sprintf('LFR: %d spikes',od.fsi_res.near_lfr_spec{iF}.spk_count), ...
            sprintf('HFR: %d spikes',od.fsi_res.near_hfr_spec{iF}.spk_count)});
        leg1.FontSize = 17;
       
        % Find PPC_peaks  
        lf = find(od.fsi_res.near_spec{iF}.freqs >= hg(1), 1, 'first');
        rf = find(od.fsi_res.near_spec{iF}.freqs <= hg(2), 1, 'last');
        % HFR Peak
        pks = findpeaks(od.fsi_res.near_hfr_spec{iF}.subsampled_ppc(lf:rf)); 
        [~, loc1] = max(od.fsi_res.near_hfr_spec{iF}.subsampled_ppc(lf+pks.loc-1));
        near_hfr_ppc_pk = pks.loc(loc1) + hg(1) - 1;
        
        % LFR Peak
        pks = findpeaks(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc(lf:rf));
        max_val = max(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc);
        [~, loc1] = max(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc(lf+pks.loc-1));
        near_lfr_ppc_pk = pks.loc(loc1) + hg(1) - 1;
        
        % Calculate mean mfr and range
        nz_trials = od.fsi_res.near_spec{iF}.mfr > 0;
        lfr_trials = od.fsi_res.near_spec{iF}.mfr <= od.fsi_res.near_spec{iF}.fr_thresh;
        lfr_trials = lfr_trials & nz_trials;
        hfr_trials = ~lfr_trials;
        hfr_trials = hfr_trials & nz_trials;
        near_lfr_spec_mean = sum(od.fsi_res.near_spec{iF}.trial_spk_count(lfr_trials))/...
            sum(od.fsi_res.near_spec{iF}.trial_spk_count(lfr_trials)'./...
            od.fsi_res.near_spec{iF}.mfr(lfr_trials));
        near_hfr_spec_mean = sum(od.fsi_res.near_spec{iF}.trial_spk_count(hfr_trials))/...
            sum(od.fsi_res.near_spec{iF}.trial_spk_count(hfr_trials)'./...
            od.fsi_res.near_spec{iF}.mfr(hfr_trials));  
          
        % Plot PPC
        ax2 = subplot(3,1,2);
        hold on;
        p3 = plot(ax2, od.fsi_res.near_spec{iF}.freqs, ...
            od.fsi_res.near_lfr_spec{iF}.subsampled_ppc, 'Color', 'red');
        q0 = xline(near_lfr_ppc_pk, 'Color', 'red', 'LineStyle','--');
        q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        p4 = plot(ax2, od.fsi_res.near_spec{iF}.freqs, ...
            od.fsi_res.near_hfr_spec{iF}.subsampled_ppc, 'Color', 'green');
        q0 = xline(near_hfr_ppc_pk, 'Color', 'green', 'LineStyle', '--');
        q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        ax2.Box = 'off';
        ax2.YTick = [];
        ax2.Title.String = 'Pairwise Phase Consistency (PPC)';
        ax2.XLabel.String = 'Frequency (in Hz)';
        ax2.Title.FontSize = 20;
        ax2.XAxis.FontSize = 18;
        ax2.XLim = hg;
        ax2.XTick = [5 10 20 30 40 50 60 70 80 90 100];
        leg2 = legend({sprintf('LFR Peak Frequency: %.2f Hz\nMean Firing Rate: %.2f Hz',...
            near_lfr_ppc_pk, near_lfr_spec_mean), ...
           sprintf('HFR Peak Frequency: %.2f Hz\nMean Firing RateL %.2f Hz',...
           near_hfr_ppc_pk, near_hfr_spec_mean)});
        leg2.FontSize = 12;
        
        %Plot PPC dif (as sanity check)
        ax3 = subplot(3,1,3);
        hold on;
        p1_p2_ppc = od.fsi_res.near_p1_spec{iF}.subsampled_ppc - ...
             od.fsi_res.near_p2_spec{iF}.subsampled_ppc;
        mean_ppc = mean(p1_p2_ppc, 1);
        sd_ppc = std(p1_p2_ppc, 1);  
        hfr_lfr_ppc = od.fsi_res.near_hfr_spec{iF}.subsampled_ppc - ...
             od.fsi_res.near_lfr_spec{iF}.subsampled_ppc;
        p5 = plot(od.fsi_res.near_spec{iF}.freqs, mean_ppc + 3*sd_ppc, '-g');
        p6 = plot(od.fsi_res.near_spec{iF}.freqs, mean_ppc - 3*sd_ppc, '-r');
        p7 = plot(od.fsi_res.near_spec{iF}.freqs, hfr_lfr_ppc, '--blue');
        p5.Annotation.LegendInformation.IconDisplayStyle = 'off';
        p6.Annotation.LegendInformation.IconDisplayStyle = 'off';
        p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if sig_fsi_type(iC) == 1
            p8 = xline(sig_fsi_min_difs(iC), '-blue');
            leg3 = legend({sprintf('LFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                sig_fsi_neg_lfr_peak_ptile(iC)*100, sig_fsi_neg_lfr_peak_ratio(iC), sig_fsi_min_difs(iC))});
        else %sig_fsi_type == 2
            p8 = xline(sig_fsi_max_difs(iC), '-blue');
            leg3 = legend({sprintf('HFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                sig_fsi_pos_hfr_peak_ptile(iC)*100, sig_fsi_pos_hfr_peak_ratio(iC), sig_fsi_max_difs(iC))});
        end
        ax3.XLabel.String = 'Frequency (in Hz)';
        ax3.Title.FontSize = 20;
        ax3.XAxis.FontSize = 18;
        ax3.XLim = hg;
        ax3.XTick = [5 10 20 30 40 50 60 70 80 90 100];
        saveas(fig, cat(2,this_label,'_FSI.fig'));
      end
   end
end
%%
for iC = 1:length(sig_msn_labels)
   %generate correct file to read
   this_label = sig_msn_labels{iC};
   fname = split(sig_msn_labels{iC},'-');
   fname = join(fname(1:4),'-');
   fname = cat(2, fname{1}, '_ft_spec.mat');
   load(fname);
   msn_labels  = od.label(od.cell_type == 1);
   for iM = 1:length(msn_labels)
      if matches(msn_labels(iM), this_label)
        close all;
        fig = figure('WindowState', 'maximized');
        % Plot STA
        ax1 = subplot(3,1,1);
        hold on;
        p1 = plot(ax1, od.msn_res.near_spec{iM}.sta_time, ...
            od.msn_res.near_lfr_spec{iM}.sta_vals, 'Color', 'red');
        p2 = plot(ax1, od.msn_res.near_spec{iM}.sta_time, ...
            od.msn_res.near_hfr_spec{iM}.sta_vals, 'Color', 'green');
        ax1.Box = 'off';
        ax1.YTick = [];
        ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
        ax1.Title.String = 'Spike Triggered Average (STA)';
        ax1.XLabel.String = 'Time (in sec)';
        ax1.Title.FontSize = 20;
        ax1.XAxis.FontSize = 18;
        leg1 = legend({sprintf('LFR: %d spikes',od.msn_res.near_lfr_spec{iM}.spk_count), ...
            sprintf('HFR: %d spikes',od.msn_res.near_hfr_spec{iM}.spk_count)});
        leg1.FontSize = 17;
       
        % Find PPC_peaks  
        lf = find(od.msn_res.near_spec{iM}.freqs >= hg(1), 1, 'first');
        rf = find(od.msn_res.near_spec{iM}.freqs <= hg(2), 1, 'last');
        % HFR Peak
        pks = findpeaks(od.msn_res.near_hfr_spec{iM}.ppc(lf:rf)); 
        [~, loc1] = max(od.msn_res.near_hfr_spec{iM}.ppc(lf+pks.loc-1));
        near_hfr_ppc_pk = pks.loc(loc1) + hg(1) - 1;
        
        % LFR Peak
        pks = findpeaks(od.msn_res.near_lfr_spec{iM}.ppc(lf:rf));
        max_val = max(od.msn_res.near_lfr_spec{iM}.ppc);
        [~, loc1] = max(od.msn_res.near_lfr_spec{iM}.ppc(lf+pks.loc-1));
        near_lfr_ppc_pk = pks.loc(loc1) + hg(1) - 1;
        
        % Calculate mean mfr and range
        nz_trials = od.msn_res.near_spec{iM}.mfr > 0;
        lfr_trials = od.msn_res.near_spec{iM}.mfr <= od.msn_res.near_spec{iM}.fr_thresh;
        lfr_trials = lfr_trials & nz_trials;
        hfr_trials = ~lfr_trials;
        hfr_trials = hfr_trials & nz_trials;
        near_lfr_spec_mean = sum(od.msn_res.near_spec{iM}.trial_spk_count(lfr_trials))/...
            sum(od.msn_res.near_spec{iM}.trial_spk_count(lfr_trials)'./...
            od.msn_res.near_spec{iM}.mfr(lfr_trials));
        near_hfr_spec_mean = sum(od.msn_res.near_spec{iM}.trial_spk_count(hfr_trials))/...
            sum(od.msn_res.near_spec{iM}.trial_spk_count(hfr_trials)'./...
            od.msn_res.near_spec{iM}.mfr(hfr_trials)); 
        
        % Plot PPC
        ax2 = subplot(3,1,2);
        hold on;
        p3 = plot(ax2, od.msn_res.near_spec{iM}.freqs, ...
            od.msn_res.near_lfr_spec{iM}.ppc, 'Color', 'red');
        q0 = xline(near_lfr_ppc_pk, 'Color', 'red', 'LineStyle','--');
        q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        p4 = plot(ax2, od.msn_res.near_spec{iM}.freqs, ...
            od.msn_res.near_hfr_spec{iM}.ppc, 'Color', 'green');
        q0 = xline(near_hfr_ppc_pk, 'Color', 'green', 'LineStyle', '--');
        q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
        ax2.Box = 'off';
        ax2.YTick = [];
        ax2.Title.String = 'Pairwise Phase Consistency (PPC)';
        ax2.XLabel.String = 'Frequency (in Hz)';
        ax2.Title.FontSize = 20;
        ax2.XAxis.FontSize = 18;
        ax2.XLim = hg;
        ax2.XTick = [5 10 20 30 40 50 60 70 80 90 100];
        leg2 = legend({sprintf('LFR Peak Frequency: %.2f Hz\nMean Firing Rate: %.2f Hz',...
            near_lfr_ppc_pk, near_lfr_spec_mean), ...
           sprintf('HFR Peak Frequency: %.2f Hz\nMean Firing RateL %.2f Hz',...
           near_hfr_ppc_pk, near_hfr_spec_mean)});
        leg2.FontSize = 12;
               
        %Plot PPC dif (as sanity check)
        ax3 = subplot(3,1,3);
        hold on;
        p1_p2_ppc = od.msn_res.near_p1_spec{iM}.ppc - ...
             od.msn_res.near_p2_spec{iM}.ppc;
        mean_ppc = mean(p1_p2_ppc, 1);
        sd_ppc = std(p1_p2_ppc, 1);  
        hfr_lfr_ppc = od.msn_res.near_hfr_spec{iM}.ppc - ...
             od.msn_res.near_lfr_spec{iM}.ppc;
        hfr_lfr_ppc = hfr_lfr_ppc';
        p5 = plot(od.msn_res.near_spec{iM}.freqs, mean_ppc + 3*sd_ppc, '-g');
        p6 = plot(od.msn_res.near_spec{iM}.freqs, mean_ppc - 3*sd_ppc, '-r');
        p7 = plot(od.msn_res.near_spec{iM}.freqs, hfr_lfr_ppc, '--blue');
        p5.Annotation.LegendInformation.IconDisplayStyle = 'off';
        p6.Annotation.LegendInformation.IconDisplayStyle = 'off';
        p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if sig_msn_type(iC) == 1
            p8 = xline(sig_msn_min_difs(iC), '-blue');
              leg3 = legend({sprintf('LFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                sig_msn_neg_lfr_peak_ptile(iC)*100, sig_msn_neg_lfr_peak_ratio(iC), sig_msn_min_difs(iC))});
        elseif sig_msn_type(iC) == 2
            p8 = xline(sig_msn_max_difs(iC), '-blue');
            leg3 = legend({sprintf('HFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                sig_msn_pos_hfr_peak_ptile(iC)*100, sig_msn_pos_hfr_peak_ratio(iC), sig_msn_max_difs(iC))});
        else %sig_msn_type == 3
            p8 = xline(sig_msn_max_difs(iC), '-blue');
            p9 = xline(sig_msn_min_difs(iC), '-black');
            leg3 = legend({sprintf('HFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                sig_msn_pos_hfr_peak_ptile(iC)*100, sig_msn_pos_hfr_peak_ratio(iC), sig_msn_max_difs(iC)),...
                sprintf('LFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                sig_msn_neg_lfr_peak_ptile(iC)*100, sig_msn_neg_lfr_peak_ratio(iC), sig_msn_min_difs(iC))});
        end        
        ax3.XLabel.String = 'Frequency (in Hz)';
        ax3.Title.FontSize = 20;
        ax3.XAxis.FontSize = 18;
        ax3.XLim = hg;
        ax3.XTick = [5 10 20 30 40 50 60 70 80 90 100];
        saveas(fig, cat(2,this_label,'_MSN.fig'));
      end
   end
end


