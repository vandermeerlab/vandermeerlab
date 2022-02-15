cd('D:\RandomVstrAnalysis\ft_results');
% cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_results');

rats = {'R117','R119','R131','R132'};
hg = [5,100];
pk_thresh = -1;
num_control_splits = 100;

load('./cellsOfInterest_cw.mat');

%Generate FSI_stuff
for iC = 1:length(fsi_cw)
   %generate correct file to read
   this_label = cat(2, fsi_cw{iC}, '.t');
   fname = split(fsi_cw{iC},'-');
   fname = join(fname(1:4),'-');
   fname = cat(2, fname{1}, '_ft_spec.mat');
   load(fname);
   fsi_labels  = od.label(od.cell_type == 2);
   for iF = 1:length(fsi_labels)
      if matches(fsi_labels(iF), this_label)
        close all;
        fig = figure('WindowState', 'maximized');
        % Plot STA
        ax1 = subplot(2,1,1);
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
        leg1.FontSize = 18;
       
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
        
        % Plot PPC
        ax2 = subplot(2,1,2);
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
        leg2 = legend({sprintf('LFR Peak Frequency: %.2f Hz',near_lfr_ppc_pk), ...
           sprintf('HFR Peak Frequency: %.2f Hz',near_hfr_ppc_pk)});
        leg2.FontSize = 18;
        saveas(fig, cat(2,this_label,'_FSI.svg'));
      end
   end
end


%Generate MSN_stuff
for iC = 1:length(msn_cw)
   %generate correct file to read
   this_label = cat(2, msn_cw{iC}, '.t');
   fname = split(msn_cw{iC},'-');
   fname = join(fname(1:4),'-');
   fname = cat(2, fname{1}, '_ft_spec.mat');
   load(fname);
   msn_labels  = od.label(od.cell_type == 1);
   for iM = 1:length(msn_labels)
      if matches(msn_labels(iM), this_label)
        close all;
        fig = figure('WindowState', 'maximized');
        % Plot STA
        ax1 = subplot(2,1,1);
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
        leg1.FontSize = 18;
       
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
        
        % Plot PPC
        ax2 = subplot(2,1,2);
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
        leg2 = legend({sprintf('LFR Peak Frequency: %.2f Hz',near_lfr_ppc_pk), ...
           sprintf('HFR Peak Frequency: %.2f Hz',near_hfr_ppc_pk)});
        leg2.FontSize = 18;
        saveas(fig, cat(2,this_label,'_MSN.svg'));
      end
   end
end


