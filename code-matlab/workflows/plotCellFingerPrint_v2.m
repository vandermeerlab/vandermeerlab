%% Script to generate cell_finger prints AND summary 
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
odir = 'D:\RandomVstrAnalysis\CellFingerPrint\';
clean_msn = 0;
clean_fsi = 0;
c1 = [75/255 0/255 146/255];  % Violet/Purple
c2 = [26/255 255/255 26/255]; % Green
c3 = [0.7 0.7 0.7]; % Gray
min_freq = 2; % Minimum frequency in Hertz
p_thresh = 99; % Percentile threshold to establish significance of phase locking
headers = {'label','lfr_min','lfr_max', 'lfr_mean', 'hfr_min', 'hfr_max', ...
    'hfr_mean', 'lfr_sts_peak','hfr_sts_peak', 'lfr_ppc_peak', 'hfr_ppc_peak'};
msn_summary = cell2table(cell(0,length(headers)), 'VariableNames', headers);
fsi_summary = cell2table(cell(0,length(headers)), 'VariableNames', headers);
%%
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        %do fsi stuff
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ~od.fsi_res.near_spec{iC}.flag_no_control_split
                %do fsi_stuff
                clean_fsi = clean_fsi + 1;
                x1 = find(round(od.fsi_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_freqs = od.fsi_res.near_spec{iC}.freqs(x1:end);

                nz_trials = od.fsi_res.near_spec{iC}.mfr > 0;
                lfr_trials = od.fsi_res.near_spec{iC}.mfr <= od.fsi_res.near_spec{iC}.fr_thresh;
                lfr_trials = lfr_trials & nz_trials;
                hfr_trials = ~lfr_trials;
                hfr_trials = hfr_trials & nz_trials;
                this_lfr_min = min(od.fsi_res.near_spec{iC}.mfr(lfr_trials));
                this_lfr_max = max(od.fsi_res.near_spec{iC}.mfr(lfr_trials));
                this_lfr_mean = sum(od.fsi_res.near_spec{iC}.trial_spk_count(lfr_trials))/ ...
                    sum(od.fsi_res.near_spec{iC}.trial_spk_count(lfr_trials)'./od.fsi_res.near_spec{iC}.mfr(lfr_trials));
                this_hfr_min = min(od.fsi_res.near_spec{iC}.mfr(hfr_trials));
                this_hfr_max = max(od.fsi_res.near_spec{iC}.mfr(hfr_trials));
                this_hfr_mean =  sum(od.fsi_res.near_spec{iC}.trial_spk_count(hfr_trials))/ ...
                    sum(od.fsi_res.near_spec{iC}.trial_spk_count(hfr_trials)'./od.fsi_res.near_spec{iC}.mfr(hfr_trials));
                
                % Plot STA
                fig = figure('WindowState', 'maximized');
                ax1 = subplot(1,3,1);
                hold on;
                p1 = plot(ax1, od.fsi_res.near_spec{iC}.sta_time, ...
                    od.fsi_res.near_lfr_spec{iC}.sta_vals, 'Color', c1);
                p2 = plot(ax1, od.fsi_res.near_spec{iC}.sta_time, ...
                    od.fsi_res.near_hfr_spec{iC}.sta_vals, 'Color', c2);
                p3 = plot(ax1, od.fsi_res.near_spec{iC}.sta_time, ...
                    od.fsi_res.near_spec{iC}.sta_vals, 'Color', c3);
%                 p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
                       
                ax1.Box = 'off';
                ax1.YTick = [];
                ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
                ax1.Title.String = 'STA';
                ax1.XLabel.String = 'Time (sec)';
                ax1.Title.FontSize = 25;
                ax1.XAxis.FontSize = 18;
                ax1.YAxis.FontSize = 18;
                ax1.XAxis.FontWeight = 'normal';
                ax1.YAxis.FontWeight = 'normal';
                ax1.TickDir = 'out';
                leg1 = legend({sprintf('LFR trials: %d spikes',od.fsi_res.near_lfr_spec{iC}.spk_count), ...
                    sprintf('HFR trials: %d spikes',od.fsi_res.near_hfr_spec{iC}.spk_count), ...
                    sprintf('All Trials: %d spikes',od.fsi_res.near_spec{iC}.spk_count)}, 'Location','best');
%                 leg1.FontSize = 17;

                % Plot STS
                ax2 = subplot(1,3,2);
                hold on;
                p4 = plot(ax2, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_lfr_spec{iC}.subsampled_sts, 'Color', c1);
                p5 = plot(ax2, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_hfr_spec{iC}.subsampled_sts, 'Color', c2);
                p6 = plot(ax2, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_spec{iC}.subsampled_sts, 'Color', c3);

                % Plot the STS thresholds
                pct_sts = prctile(od.fsi_res.near_spec{iC}.shuf_sts, p_thresh);
                plot(ax2, od.fsi_res.near_spec{iC}.freqs, pct_sts, '--black');
                lfr_sts_mask = od.fsi_res.near_lfr_spec{iC}.subsampled_sts>pct_sts;
                hfr_sts_mask = od.fsi_res.near_hfr_spec{iC}.subsampled_sts>pct_sts;

                [~, this_lfr_sts_max] = max(od.fsi_res.near_lfr_spec{iC}.subsampled_sts(x1:end).*(lfr_sts_mask(x1:end)));
                if ~isempty(this_lfr_sts_max) 
                    this_lfr_sts_max = this_freqs(this_lfr_sts_max);
                    xline(ax2, this_lfr_sts_max, 'Color', c1,  'LineStyle', '--');
                else
                    this_lfr_sts_max = NaN;
                end
                [~, this_hfr_sts_max] = max(od.fsi_res.near_hfr_spec{iC}.subsampled_sts(x1:end).*(hfr_sts_mask(x1:end)));
                if ~isempty(this_hfr_sts_max) 
                    this_hfr_sts_max = this_freqs(this_hfr_sts_max);
                    xline(ax2, this_hfr_sts_max, 'Color', c2,  'LineStyle', '--');
                else
                    this_hfr_sts_max = NaN;
                end

                ax2.Box = 'off';
                ax2.YTick = [];
                ax2.XLim = [x1 100];
                ax2.XTick = [x1 10 25 50 75 100];
                ax2.Title.String = 'STS';
                ax2.XLabel.String = 'Frequency (Hz)';
                ax2.Title.FontSize = 25;
                ax2.XAxis.FontSize = 18;
                ax2.XAxis.FontWeight = 'normal';
                ax2.TickDir = 'out';
                leg2 = legend({'LFR','HFR','All Trials', sprintf('%2d thresh', p_thresh)}, 'Location','best');
                
                % Plot PPC
                ax3 = subplot(1,3,3);
                hold on;
                p9 = plot(ax3, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_lfr_spec{iC}.subsampled_ppc, 'Color', c1);
                p10 = plot(ax3, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_hfr_spec{iC}.subsampled_ppc, 'Color', c2);
                p11 = plot(ax3, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_spec{iC}.subsampled_ppc, 'Color', c3);

                % Plot the PPC thresholds
                pct_ppc = prctile(od.fsi_res.near_spec{iC}.shuf_ppc, p_thresh);
                plot(ax3, od.fsi_res.near_spec{iC}.freqs, pct_ppc, '--black');
                lfr_ppc_mask = od.fsi_res.near_lfr_spec{iC}.subsampled_ppc>pct_ppc;
                hfr_ppc_mask = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc>pct_ppc;

                [~, this_lfr_ppc_max] = max(od.fsi_res.near_lfr_spec{iC}.subsampled_ppc(x1:end).*(lfr_ppc_mask(x1:end)));
                if ~isempty(this_lfr_ppc_max) 
                    this_lfr_ppc_max = this_freqs(this_lfr_ppc_max);
                    xline(ax3, this_lfr_ppc_max, 'Color', c1,  'LineStyle', '--');
                else
                    this_lfr_ppc_max = NaN;
                end
                [~, this_hfr_ppc_max] = max(od.fsi_res.near_hfr_spec{iC}.subsampled_ppc(x1:end).*(hfr_ppc_mask(x1:end)));
                if ~isempty(this_hfr_ppc_max) 
                    this_hfr_ppc_max = this_freqs(this_hfr_ppc_max);
                    xline(ax3, this_hfr_ppc_max, 'Color', c2,  'LineStyle', '--');
                else
                    this_hfr_ppc_max = NaN;
                end

                ax3.Box = 'off';
                ax3.XLim = [x1 100];
                ax3.XTick = [x1 10 25 50 75 100];
                ax3.Title.String = 'PPC';
                ax3.XLabel.String = 'Frequency (Hz)';
                ax3.Title.FontSize = 25;
                ax3.XAxis.FontSize = 18;
                ax3.XAxis.FontWeight = 'normal';
                ax3.YAxis.FontSize = 18;
                ax3.YAxis.FontWeight = 'normal';
                ax3.TickDir = 'out';
                leg3 = legend({sprintf('%.2f Hz', this_lfr_mean), sprintf('%.2f Hz', this_hfr_mean), ...
                    'All Trials', sprintf('%2d thresh', p_thresh)}, 'Location','best');
                
                this_row = {fsi_labels{iC}, this_lfr_min, this_lfr_max, ...
                    this_lfr_mean, this_hfr_min, this_hfr_max, this_hfr_mean, ...
                    this_lfr_sts_max, this_hfr_sts_max, this_lfr_ppc_max, this_hfr_ppc_max};
                
                fsi_summary = [fsi_summary; this_row];
	            print(fig, '-dpng', '-r300', strcat(odir, fsi_labels{iC}, '_FSI'));
                close;
            end
        end

        %do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                %do msn_stuff
                clean_msn = clean_msn + 1;
                                x1 = find(round(od.msn_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                this_freqs = od.msn_res.near_spec{iC}.freqs(x1:end);

                nz_trials = od.msn_res.near_spec{iC}.mfr > 0;
                lfr_trials = od.msn_res.near_spec{iC}.mfr <= od.msn_res.near_spec{iC}.fr_thresh;
                lfr_trials = lfr_trials & nz_trials;
                hfr_trials = ~lfr_trials;
                hfr_trials = hfr_trials & nz_trials;
                this_lfr_min = min(od.msn_res.near_spec{iC}.mfr(lfr_trials));
                this_lfr_max = max(od.msn_res.near_spec{iC}.mfr(lfr_trials));
                this_lfr_mean = sum(od.msn_res.near_spec{iC}.trial_spk_count(lfr_trials))/ ...
                    sum(od.msn_res.near_spec{iC}.trial_spk_count(lfr_trials)'./od.msn_res.near_spec{iC}.mfr(lfr_trials));
                this_hfr_min = min(od.msn_res.near_spec{iC}.mfr(hfr_trials));
                this_hfr_max = max(od.msn_res.near_spec{iC}.mfr(hfr_trials));
                this_hfr_mean =  sum(od.msn_res.near_spec{iC}.trial_spk_count(hfr_trials))/ ...
                    sum(od.msn_res.near_spec{iC}.trial_spk_count(hfr_trials)'./od.msn_res.near_spec{iC}.mfr(hfr_trials));
                
                % Plot STA
                fig = figure('WindowState', 'maximized');
                ax1 = subplot(1,3,1);
                hold on;
                p1 = plot(ax1, od.msn_res.near_spec{iC}.sta_time, ...
                    od.msn_res.near_lfr_spec{iC}.sta_vals, 'Color', c1);
                p2 = plot(ax1, od.msn_res.near_spec{iC}.sta_time, ...
                    od.msn_res.near_hfr_spec{iC}.sta_vals, 'Color', c2);
                p3 = plot(ax1, od.msn_res.near_spec{iC}.sta_time, ...
                    od.msn_res.near_spec{iC}.sta_vals, 'Color', c3);
%                 p3.Annotation.LegendInformation.IconDisplayStyle = 'off';
                       
                ax1.Box = 'off';
                ax1.YTick = [];
                ax1.XTick = [-0.5 -0.25 0 0.25 0.5];
                ax1.Title.String = 'STA';
                ax1.XLabel.String = 'Time (sec)';
                ax1.Title.FontSize = 25;
                ax1.XAxis.FontSize = 18;
                ax1.YAxis.FontSize = 18;
                ax1.XAxis.FontWeight = 'normal';
                ax1.YAxis.FontWeight = 'normal';
                ax1.TickDir = 'out';
                leg1 = legend({sprintf('LFR trials: %d spikes',od.msn_res.near_lfr_spec{iC}.spk_count), ...
                    sprintf('HFR trials: %d spikes',od.msn_res.near_hfr_spec{iC}.spk_count), ...
                    sprintf('All Trials: %d spikes',od.msn_res.near_spec{iC}.spk_count)}, 'Location','best');
%                 leg1.FontSize = 17;

                % Plot STS
                ax2 = subplot(1,3,2);
                hold on;
                p4 = plot(ax2, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_lfr_spec{iC}.sts_vals, 'Color', c1);
                p5 = plot(ax2, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_hfr_spec{iC}.sts_vals, 'Color', c2);
                p6 = plot(ax2, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_spec{iC}.sts_vals, 'Color', c3);

                % Plot the STS thresholds
                pct_sts = prctile(od.msn_res.near_spec{iC}.shuf_sts, p_thresh);
                plot(ax2, od.msn_res.near_spec{iC}.freqs, pct_sts, '--black');
                lfr_sts_mask = od.msn_res.near_lfr_spec{iC}.sts_vals>pct_sts;
                hfr_sts_mask = od.msn_res.near_hfr_spec{iC}.sts_vals>pct_sts;

                [~, this_lfr_sts_max] = max(od.msn_res.near_lfr_spec{iC}.sts_vals(x1:end).*(lfr_sts_mask(x1:end)));
                if ~isempty(this_lfr_sts_max) 
                    this_lfr_sts_max = this_freqs(this_lfr_sts_max);
                    xline(ax2, this_lfr_sts_max, 'Color', c1,  'LineStyle', '--');
                else
                    this_lfr_sts_max = NaN;
                end
                [~, this_hfr_sts_max] = max(od.msn_res.near_hfr_spec{iC}.sts_vals(x1:end).*(hfr_sts_mask(x1:end)));
                if ~isempty(this_hfr_sts_max) 
                    this_hfr_sts_max = this_freqs(this_hfr_sts_max);
                    xline(ax2, this_hfr_sts_max, 'Color', c2,  'LineStyle', '--');
                else
                    this_hfr_sts_max = NaN;
                end

                ax2.Box = 'off';
                ax2.YTick = [];
                ax2.XLim = [x1 100];
                ax2.XTick = [x1 10 25 50 75 100];
                ax2.Title.String = 'STS';
                ax2.XLabel.String = 'Frequency (Hz)';
                ax2.Title.FontSize = 25;
                ax2.XAxis.FontSize = 18;
                ax2.XAxis.FontWeight = 'normal';
                ax2.TickDir = 'out';
                leg2 = legend({'LFR','HFR','All Trials', sprintf('%2d thresh', p_thresh)}, 'Location','best');
                
                % Plot PPC
                ax3 = subplot(1,3,3);
                hold on;
                p9 = plot(ax3, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_lfr_spec{iC}.ppc', 'Color', c1);
                p10 = plot(ax3, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_hfr_spec{iC}.ppc', 'Color', c2);
                p11 = plot(ax3, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_spec{iC}.ppc', 'Color', c3);

                % Plot the PPC thresholds
                pct_ppc = prctile(od.msn_res.near_spec{iC}.shuf_ppc, p_thresh);
                plot(ax3, od.msn_res.near_spec{iC}.freqs, pct_ppc, '--black');
                lfr_ppc_mask = od.msn_res.near_lfr_spec{iC}.ppc'>pct_ppc;
                hfr_ppc_mask = od.msn_res.near_hfr_spec{iC}.ppc'>pct_ppc;

                [~, this_lfr_ppc_max] = max(od.msn_res.near_lfr_spec{iC}.ppc(x1:end)'.*(lfr_ppc_mask(x1:end)));
                if ~isempty(this_lfr_ppc_max) 
                    this_lfr_ppc_max = this_freqs(this_lfr_ppc_max);
                    xline(ax3, this_lfr_ppc_max, 'Color', c1,  'LineStyle', '--');
                else
                    this_lfr_ppc_max = NaN;
                end
                [~, this_hfr_ppc_max] = max(od.msn_res.near_hfr_spec{iC}.ppc(x1:end)'.*(hfr_ppc_mask(x1:end)));
                if ~isempty(this_hfr_ppc_max) 
                    this_hfr_ppc_max = this_freqs(this_hfr_ppc_max);
                    xline(ax3, this_hfr_ppc_max, 'Color', c2,  'LineStyle', '--');
                else
                    this_hfr_ppc_max = NaN;
                end

                ax3.Box = 'off';
                ax3.XLim = [x1 100];
                ax3.XTick = [x1 10 25 50 75 100];
                ax3.Title.String = 'PPC';
                ax3.XLabel.String = 'Frequency (Hz)';
                ax3.Title.FontSize = 25;
                ax3.XAxis.FontSize = 18;
                ax3.XAxis.FontWeight = 'normal';
                ax3.YAxis.FontSize = 18;
                ax3.YAxis.FontWeight = 'normal';
                ax3.TickDir = 'out';
                leg3 = legend({sprintf('%.2f Hz', this_lfr_mean), sprintf('%.2f Hz', this_hfr_mean), ...
                    'All Trials', sprintf('%2d thresh', p_thresh)}, 'Location','best');
                
                this_row = {msn_labels{iC}, this_lfr_min, this_lfr_max, ...
                    this_lfr_mean, this_hfr_min, this_hfr_max, this_hfr_mean, ...
                    this_lfr_sts_max, this_hfr_sts_max, this_lfr_ppc_max, this_hfr_ppc_max};
                
                msn_summary = [msn_summary; this_row];
                print(fig, '-dpng', '-r300', strcat(odir, msn_labels{iC}, '_MSN'));
                close;
            end
        end
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
writetable(msn_summary, 'msn_summary.xls');
writetable(fsi_summary, 'fsi_summary.xls');