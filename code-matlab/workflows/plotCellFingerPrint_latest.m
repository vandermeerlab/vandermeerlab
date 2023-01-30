cd('E:\Dropbox (Dartmouth College)\AnalysisResults\FieldTripResults\ft_results');
rats = {'R117','R119','R131','R132'};
fw = [2,100];
pk_thresh = -1;
num_control_splits = 100;

fsi_mfr_difs = [];
fsi_freq_difs = [];
msn_mfr_difs = [];
msn_freq_difs = [];
load('cellsOfInterest_cw.mat');
load('significant_dif.mat');

set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesFontWeight','bold')
%%
%Generate FSI_stuff
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)        
        load(ofiles(jdx).name);
        fsi_labels  = od.label(od.cell_type == 2);
        for iF = 1:length(fsi_labels)
            this_label = extractBefore(fsi_labels{iF}, '.t');
            if sum(matches(fsi_counter, this_label)) == 0
                continue;
            end
            fig = figure('WindowState', 'maximized');
            % Plot STA
            ax1 = subplot(1,1,1);
            hold on;
            p1 = plot(ax1, od.fsi_res.near_spec{iF}.sta_time, ...
                od.fsi_res.near_lfr_spec{iF}.sta_vals, 'Color', 'red');
            p2 = plot(ax1, od.fsi_res.near_spec{iF}.sta_time, ...
                od.fsi_res.near_hfr_spec{iF}.sta_vals, 'Color', 'green');
            p3 = plot(ax1, od.fsi_res.near_spec{iF}.sta_time, ...
                od.fsi_res.near_spec{iF}.sta_vals, 'Color', 'blue');
%             p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
                       
            ax1.Box = 'off';
            ax1.YTick = [];
            ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
            ax1.Title.String = 'STA';
            ax1.XLabel.String = 'Time (in sec)';
            ax1.Title.FontSize = 25;
            ax1.XAxis.FontSize = 18;
            ax1.YAxis.FontSize = 18;
            ax1.XAxis.FontWeight = 'normal';
            ax1.YAxis.FontWeight = 'normal';
            ax1.TickDir = 'out';
            leg1 = legend({sprintf('LFR trials: %d spikes',od.fsi_res.near_lfr_spec{iF}.spk_count), ...
                sprintf('HFR trials: %d spikes',od.fsi_res.near_hfr_spec{iF}.spk_count), ...
                sprintf('All Trials: %d spikes',od.fsi_res.near_spec{iF}.spk_count)});
            leg1.FontSize = 17;

            % Find PPC_peaks  
            lf = find(od.fsi_res.near_spec{iF}.freqs >= fw(1), 1, 'first');
            rf = find(od.fsi_res.near_spec{iF}.freqs <= fw(2), 1, 'last');
            
            % All Trials Peak
            pks = findpeaks(od.fsi_res.near_spec{iF}.subsampled_ppc(lf:rf)); 
            [~, loc1] = max(od.fsi_res.near_spec{iF}.subsampled_ppc(lf+pks.loc-1));
            near_ppc_pk = od.fsi_res.near_spec{iF}.freqs(pks.loc(loc1) + fw(1) - 1);
            
            % HFR Peak
            pks = findpeaks(od.fsi_res.near_hfr_spec{iF}.subsampled_ppc(lf:rf)); 
            [~, loc1] = max(od.fsi_res.near_hfr_spec{iF}.subsampled_ppc(lf+pks.loc-1));
            near_hfr_ppc_pk = od.fsi_res.near_spec{iF}.freqs(pks.loc(loc1) + fw(1) - 1);

            % LFR Peak
            pks = findpeaks(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc(lf:rf));
            max_val = max(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc);
            [~, loc1] = max(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc(lf+pks.loc-1));
            near_lfr_ppc_pk = od.fsi_res.near_spec{iF}.freqs(pks.loc(loc1) + fw(1) - 1);

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
            % Add 'SIG' to filenames if this is significant
            if sig_dif_fsi(matches(fsi_counter, this_label)) ~= 0
                saveas(fig, cat(2,this_label,'_SIG_STA_FSI.fig'));
            else
                saveas(fig, cat(2,this_label,'_STA_FSI.fig'));
            end
            close all;

            % Plot PPC
            fig = figure('WindowState', 'maximized');
            ax2 = subplot(1,1,1);
            hold on;
            p4 = plot(ax2, od.fsi_res.near_spec{iF}.freqs, ...
                od.fsi_res.near_lfr_spec{iF}.subsampled_ppc, 'Color', 'red');
            q0 = xline(near_lfr_ppc_pk, 'Color', 'red', 'LineStyle','--');
            q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p5 = plot(ax2, od.fsi_res.near_spec{iF}.freqs, ...
                od.fsi_res.near_hfr_spec{iF}.subsampled_ppc, 'Color', 'green');
            q0 = xline(near_hfr_ppc_pk, 'Color', 'green', 'LineStyle', '--');
            q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p6 = plot(ax2, od.fsi_res.near_spec{iF}.freqs, ...
                od.fsi_res.near_spec{iF}.subsampled_ppc, 'Color', 'blue');
            q0 = xline(near_ppc_pk, 'Color', 'blue', 'LineStyle', '--');
            q0.Annotation.LegendInformation.IconDisplayStyle = 'off';

            ax2.Box = 'off';
            ax2.Title.String = 'PPC';
%             ax2.YLabel.String = 'PPC';
            ax2.XLabel.String = 'Frequency (in Hz)';
            ax2.Title.FontSize = 25;
            ax2.XAxis.FontSize = 18;
            ax2.XLim = fw;
%             ax2.YLim = [0 0.04];
            ax2.XTick = [5 10 20 30 40 50 60 70 80 90 100];
            ax2.XAxis.FontWeight = 'normal';
            ax2.YAxis.FontWeight = 'normal';
            ax2.TickDir = 'out';
            leg2 = legend({sprintf('LFR Peak Frequency: %.2f Hz\nMean Firing Rate: %.2f Hz',...
                near_lfr_ppc_pk, near_lfr_spec_mean), ...
               sprintf('HFR Peak Frequency: %.2f Hz\nMean Firing Rate %.2f Hz',...
               near_hfr_ppc_pk, near_hfr_spec_mean)});
            leg2.FontSize = 12;
            fsi_mfr_difs = [fsi_mfr_difs near_hfr_spec_mean-near_lfr_spec_mean];
            fsi_freq_difs = [fsi_freq_difs near_hfr_ppc_pk-near_lfr_ppc_pk];

            % Add 'SIG' to filenames if this is significant
            if sig_dif_fsi(matches(fsi_counter, this_label)) ~= 0
                saveas(fig, cat(2,this_label,'_SIG_PPC_FSI.fig'));
            else
                saveas(fig, cat(2,this_label,'_PPC_FSI.fig'));
            end
            close all;

            %Plot PPC dif (as sanity check)
            fig = figure('WindowState', 'maximized');
            ax3 = subplot(1,1,1);
            hold on;
            p1_p2_ppc = od.fsi_res.near_p1_spec{iF}.subsampled_ppc - ...
                 od.fsi_res.near_p2_spec{iF}.subsampled_ppc;
            mean_ppc = mean(p1_p2_ppc, 1);
            sd_ppc = std(p1_p2_ppc, 1);  
            hfr_lfr_ppc = od.fsi_res.near_hfr_spec{iF}.subsampled_ppc - ...
                 od.fsi_res.near_lfr_spec{iF}.subsampled_ppc;
            p7 = plot(od.fsi_res.near_spec{iF}.freqs, mean_ppc + 3*sd_ppc, '-g');
            p8 = plot(od.fsi_res.near_spec{iF}.freqs, mean_ppc - 3*sd_ppc, '-r');
            p9 = plot(od.fsi_res.near_spec{iF}.freqs, hfr_lfr_ppc, '--blue');
            p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p8.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p9.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            if sig_dif_fsi(matches(fsi_counter, this_label)) == 1
                % hacky_way to get the inde
                this_idx =  find(find(sig_dif_fsi) == find(matches(fsi_counter, this_label)));
                p8 = xline(sig_fsi_min_difs(this_idx), '-blue');
                leg3 = legend({sprintf('LFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                    sig_fsi_neg_lfr_peak_ptile(this_idx)*100, sig_fsi_neg_lfr_peak_ratio(this_idx), sig_fsi_min_difs(this_idx))});
            elseif sig_dif_fsi(matches(fsi_counter, this_label)) == 2
                this_idx =  find(find(sig_dif_fsi) == find(matches(fsi_counter, this_label)));
                p8 = xline(sig_fsi_max_difs(this_idx), '-blue');
                leg3 = legend({sprintf('HFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                    sig_fsi_pos_hfr_peak_ptile(this_idx)*100, sig_fsi_pos_hfr_peak_ratio(this_idx), sig_fsi_max_difs(this_idx))});
            end
            ax3.XLabel.String = 'Frequency (in Hz)';
            ax3.Title.FontSize = 25;
            ax3.XAxis.FontSize = 18;
            ax3.XLim = fw;
            ax3.XTick = [5 10 20 30 40 50 60 70 80 90 100];
            ax3.XAxis.FontWeight = 'normal';
            ax3.YAxis.FontWeight = 'normal';
            ax3.TickDir = 'out';
            % Add 'SIG' to filenames if this is significant
            if sig_dif_fsi(matches(fsi_counter, this_label)) ~= 0
                saveas(fig, cat(2,this_label,'_SIG_PPC_diff_FSI.fig'));
            else
                saveas(fig, cat(2,this_label,'_PPC_diff_FSI.fig'));
            end
            close all
        end
        
        msn_labels = od.label(od.cell_type == 1);      
        for iM = 1:length(msn_labels)
            this_label = extractBefore(msn_labels{iM}, '.t');
            if sum(matches(msn_counter, this_label)) == 0
                continue;
            end
            fig = figure('WindowState', 'maximized');
            % Plot STA
            ax1 = subplot(1,1,1);
            hold on;
            p1 = plot(ax1, od.msn_res.near_spec{iM}.sta_time, ...
                od.msn_res.near_lfr_spec{iM}.sta_vals, 'Color', 'red');
            p2 = plot(ax1, od.msn_res.near_spec{iM}.sta_time, ...
                od.msn_res.near_hfr_spec{iM}.sta_vals, 'Color', 'green');
            p3 = plot(ax1, od.msn_res.near_spec{iM}.sta_time, ...
                od.msn_res.near_spec{iM}.sta_vals, 'Color', 'blue');
        %             p3.Annotation.LegendInformation.IconDisplayStyle = 'off';

            ax1.Box = 'off';
            ax1.YTick = [];
            ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
            ax1.Title.String = 'STA';
            ax1.XLabel.String = 'Time (in sec)';
            ax1.Title.FontSize = 25;
            ax1.XAxis.FontSize = 18;
            ax1.YAxis.FontSize = 18;
            ax1.XAxis.FontWeight = 'normal';
            ax1.YAxis.FontWeight = 'normal';
            ax1.TickDir = 'out';
            leg1 = legend({sprintf('LFR trials: %d spikes',od.msn_res.near_lfr_spec{iM}.spk_count), ...
                sprintf('HFR trials: %d spikes',od.msn_res.near_hfr_spec{iM}.spk_count), ...
                sprintf('All Trials: %d spikes',od.msn_res.near_spec{iM}.spk_count)});
            leg1.FontSize = 17;

            % Find PPC_peaks  
            lf = find(od.msn_res.near_spec{iM}.freqs >= fw(1), 1, 'first');
            rf = find(od.msn_res.near_spec{iM}.freqs <= fw(2), 1, 'last');

            % All Trials Peak
            pks = findpeaks(od.msn_res.near_spec{iM}.ppc(lf:rf)); 
            [~, loc1] = max(od.msn_res.near_spec{iM}.ppc(lf+pks.loc-1));
            near_ppc_pk = od.msn_res.near_spec{iM}.freqs(pks.loc(loc1) + fw(1) - 1);

            % HFR Peak
            pks = findpeaks(od.msn_res.near_hfr_spec{iM}.ppc(lf:rf)); 
            [~, loc1] = max(od.msn_res.near_hfr_spec{iM}.ppc(lf+pks.loc-1));
            near_hfr_ppc_pk = od.msn_res.near_spec{iM}.freqs(pks.loc(loc1) + fw(1) - 1);

            % LFR Peak
            pks = findpeaks(od.msn_res.near_lfr_spec{iM}.ppc(lf:rf));
            max_val = max(od.msn_res.near_lfr_spec{iM}.ppc);
            [~, loc1] = max(od.msn_res.near_lfr_spec{iM}.ppc(lf+pks.loc-1));
            near_lfr_ppc_pk = od.msn_res.near_spec{iM}.freqs(pks.loc(loc1) + fw(1) - 1);

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
            % Add 'SIG' to filenames if this is significant
            if sig_dif_msn(matches(msn_counter, this_label)) ~= 0
                saveas(fig, cat(2,this_label,'_SIG_STA_MSN.fig'));
            else
                saveas(fig, cat(2,this_label,'_STA_MSN.fig'));
            end
            close all;

            % Plot PPC
            fig = figure('WindowState', 'maximized');
            ax2 = subplot(1,1,1);
            hold on;
            p4 = plot(ax2, od.msn_res.near_spec{iM}.freqs, ...
                od.msn_res.near_lfr_spec{iM}.ppc, 'Color', 'red');
            q0 = xline(near_lfr_ppc_pk, 'Color', 'red', 'LineStyle','--');
            q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p5 = plot(ax2, od.msn_res.near_spec{iM}.freqs, ...
                od.msn_res.near_hfr_spec{iM}.ppc, 'Color', 'green');
            q0 = xline(near_hfr_ppc_pk, 'Color', 'green', 'LineStyle', '--');
            q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p6 = plot(ax2, od.msn_res.near_spec{iM}.freqs, ...
                od.msn_res.near_spec{iM}.ppc, 'Color', 'blue');
            q0 = xline(near_ppc_pk, 'Color', 'blue', 'LineStyle', '--');
            q0.Annotation.LegendInformation.IconDisplayStyle = 'off';

            ax2.Box = 'off';
            ax2.Title.String = 'PPC';
        %             ax2.YLabel.String = 'PPC';
            ax2.XLabel.String = 'Frequency (in Hz)';
            ax2.Title.FontSize = 25;
            ax2.XAxis.FontSize = 18;
            ax2.XLim = fw;
        %             ax2.YLim = [0 0.04];
            ax2.XTick = [5 10 20 30 40 50 60 70 80 90 100];
            ax2.XAxis.FontWeight = 'normal';
            ax2.YAxis.FontWeight = 'normal';
            ax2.TickDir = 'out';
            leg2 = legend({sprintf('LFR Peak Frequency: %.2f Hz\nMean Firing Rate: %.2f Hz',...
                near_lfr_ppc_pk, near_lfr_spec_mean), ...
               sprintf('HFR Peak Frequency: %.2f Hz\nMean Firing Rate %.2f Hz',...
               near_hfr_ppc_pk, near_hfr_spec_mean)});
            leg2.FontSize = 12;
            msn_mfr_difs = [msn_mfr_difs near_hfr_spec_mean-near_lfr_spec_mean];
            msn_freq_difs = [msn_freq_difs near_hfr_ppc_pk-near_lfr_ppc_pk];

            % Add 'SIG' to filenames if this is significant
            if sig_dif_msn(matches(msn_counter, this_label)) ~= 0
                saveas(fig, cat(2,this_label,'_SIG_PPC_MSN.fig'));
            else
                saveas(fig, cat(2,this_label,'_PPC_MSN.fig'));
            end
            close all;

            %Plot PPC dif (as sanity check)
            fig = figure('WindowState', 'maximized');
            ax3 = subplot(1,1,1);
            hold on;
            p1_p2_ppc = od.msn_res.near_p1_spec{iM}.ppc - ...
                 od.msn_res.near_p2_spec{iM}.ppc;
            mean_ppc = mean(p1_p2_ppc, 1);
            sd_ppc = std(p1_p2_ppc, 1);  
            hfr_lfr_ppc = od.msn_res.near_hfr_spec{iM}.ppc - ...
                 od.msn_res.near_lfr_spec{iM}.ppc;
            p7 = plot(od.msn_res.near_spec{iM}.freqs, mean_ppc + 3*sd_ppc, '-g');
            p8 = plot(od.msn_res.near_spec{iM}.freqs, mean_ppc - 3*sd_ppc, '-r');
            p9 = plot(od.msn_res.near_spec{iM}.freqs, hfr_lfr_ppc, '--blue');
            p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p8.Annotation.LegendInformation.IconDisplayStyle = 'off';
            p9.Annotation.LegendInformation.IconDisplayStyle = 'off';

            if sig_dif_msn(matches(msn_counter, this_label)) == 1
                % hacky_way to get the inde
                this_idx =  find(find(sig_dif_msn) == find(matches(msn_counter, this_label)));
                p8 = xline(sig_msn_min_difs(this_idx), '-blue');
                leg3 = legend({sprintf('LFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                    sig_msn_neg_lfr_peak_ptile(this_idx)*100, sig_msn_neg_lfr_peak_ratio(this_idx), sig_msn_min_difs(this_idx))});
            elseif sig_dif_msn(matches(msn_counter, this_label)) == 2
                this_idx =  find(find(sig_dif_msn) == find(matches(msn_counter, this_label)));
                p8 = xline(sig_msn_max_difs(this_idx), '-blue');
                leg3 = legend({sprintf('HFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                    sig_msn_pos_hfr_peak_ptile(this_idx)*100, sig_msn_pos_hfr_peak_ratio(this_idx), sig_msn_max_difs(this_idx))});
            elseif sig_dif_msn(matches(msn_counter, this_label)) == 3
                this_idx =  find(find(sig_dif_msn) == find(matches(msn_counter, this_label)));
                p8 = xline(sig_msn_min_difs(this_idx), '-blue');
                p9 = xline(sig_msn_min_difs(this_idx), '-black');
                leg3 = legend({sprintf('LFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                    sig_msn_neg_lfr_peak_ptile(this_idx)*100, sig_msn_neg_lfr_peak_ratio(this_idx), sig_msn_min_difs(this_idx)), ...
                    sprintf('HFR Peak percentile: %.2f \nPeak Ratio: %.2f \nFrequency: %.2f Hz', ...
                    sig_msn_pos_hfr_peak_ptile(this_idx)*100, sig_msn_pos_hfr_peak_ratio(this_idx), sig_msn_max_difs(this_idx))});
            end
            ax3.XLabel.String = 'Frequency (in Hz)';
            ax3.Title.FontSize = 25;
            ax3.XAxis.FontSize = 18;
            ax3.XLim = fw;
            ax3.XTick = [5 10 20 30 40 50 60 70 80 90 100];
            ax3.XAxis.FontWeight = 'normal';
            ax3.YAxis.FontWeight = 'normal';
            ax3.TickDir = 'out';
            % Add 'SIG' to filenames if this is significant
            if sig_dif_msn(matches(msn_counter, this_label)) ~= 0
                saveas(fig, cat(2,this_label,'_SIG_PPC_diff_MSN.fig'));
            else
                saveas(fig, cat(2,this_label,'_PPC_diff_MSN.fig'));
            end
            close all
        end
    end
end


