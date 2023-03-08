%% Script to generate cell_finger prints
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
c1 = [75/255 0/255 146/255];  % Violet/Purple
c2 = [26/255 255/255 26/255]; % Green
c3 = [0.7 0.7 0.7]; %Gray
min_freq = 2;
% Load Significance info
load('./significance.mat');
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
                this_sig = strcmp(sig_fsi(:,1),fsi_labels{iC});
                this_sig = sig_fsi(this_sig,:);
                x1 = find(round(od.fsi_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                
                % Plot STA
                fig = figure('WindowState', 'maximized');
                ax1 = subplot(2,3,1);
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
                ax2 = subplot(2,3,2);
                hold on;
                p4 = plot(ax2, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_lfr_spec{iC}.subsampled_sts, 'Color', c1);
                p5 = plot(ax2, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_hfr_spec{iC}.subsampled_sts, 'Color', c2);
                p6 = plot(ax2, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_spec{iC}.subsampled_sts, 'Color', c3);
                % show vertical lines at points of significance
                if ~this_sig{4} == 0 
                    p7 = xline(ax2, this_sig{4}, 'LineStyle', '--', 'Color',c2);
                    p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{5} == 0
                    p8 = xline(ax2, this_sig{5}, 'LineStyle', '--', 'Color',c1);
                    p8.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
                leg2 = legend({'LFR','HFR','All Trials'}, 'Location','best');
                
                % Plot PPC
                ax3 = subplot(2,3,3);
                hold on;
                p9 = plot(ax3, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_lfr_spec{iC}.subsampled_ppc, 'Color', c1);
                p10 = plot(ax3, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_hfr_spec{iC}.subsampled_ppc, 'Color', c2);
                p11 = plot(ax3, od.fsi_res.near_spec{iC}.freqs, ...
                    od.fsi_res.near_spec{iC}.subsampled_ppc, 'Color', c3);
                % show vertical lines at points of significance
                if ~this_sig{2} == 0 
                    p12 = xline(ax3, this_sig{2}, 'LineStyle', '--', 'Color',c2);
                    p12.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{3} == 0
                    p13 = xline(ax3, this_sig{3}, 'LineStyle', '--', 'Color',c1);
                    p13.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
                leg3 = legend({'LFR','HFR','All Trials'}, 'Location','best');

                % Plot STS dif
                ax4 = subplot(2,3,5);
                hold on;
                p1_p2_sts = od.fsi_res.near_p1_spec{iC}.subsampled_sts - od.fsi_res.near_p2_spec{iC}.subsampled_sts;
                mean_sts = mean(p1_p2_sts,1);
                sd_sts = std(p1_p2_sts);
                hfr_lfr_sts = od.fsi_res.near_hfr_spec{iC}.subsampled_sts - od.fsi_res.near_lfr_spec{iC}.subsampled_sts;
                p14 = plot(ax4, od.fsi_res.near_spec{iC}.freqs, mean_sts + (3*sd_sts), '--m');
                p15 = plot(ax4, od.fsi_res.near_spec{iC}.freqs, mean_sts - (3*sd_sts), '--m');
                p15.Annotation.LegendInformation.IconDisplayStyle = 'off';
                p16 = plot(ax4, od.fsi_res.near_spec{iC}.freqs, hfr_lfr_sts, 'black');
                % show vertical lines at points of significance
                if ~this_sig{4} == 0 
                    p17 = xline(ax4, this_sig{4}, 'LineStyle', '--', 'Color',c2);
                    p17.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{5} == 0
                    p18 = xline(ax4, this_sig{5}, 'LineStyle', '--', 'Color',c1);
                    p18.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end 
                ax4.Box = 'off';
                ax4.YTick = [];
                ax4.XLim = [x1 100];
                ax4.XTick = [x1 10 25 50 75 100];
                ax4.Title.String = 'STS diff';
                ax4.XLabel.String = 'Frequency (Hz)';
                ax4.Title.FontSize = 25;
                ax4.XAxis.FontSize = 18;
                ax4.XAxis.FontWeight = 'normal';
                ax4.TickDir = 'out';
                leg4 = legend({'Mean +/- 3*SD','HFR - LFR'}, 'Location','best');


                % Plot PPC dif
                ax5 = subplot(2,3,6);
                hold on;
                p1_p2_ppc = od.fsi_res.near_p1_spec{iC}.subsampled_ppc - od.fsi_res.near_p2_spec{iC}.subsampled_ppc;
                mean_ppc = mean(p1_p2_ppc, 1);
                sd_ppc = std(p1_p2_ppc);
                hfr_lfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc - od.fsi_res.near_lfr_spec{iC}.subsampled_ppc;
                p19 = plot(ax5, od.fsi_res.near_spec{iC}.freqs, mean_ppc + (3*sd_ppc), '--m');
                p20 = plot(ax5, od.fsi_res.near_spec{iC}.freqs, mean_ppc - (3*sd_ppc), '--m');
                p20.Annotation.LegendInformation.IconDisplayStyle = 'off';
                p21 = plot(ax5, od.fsi_res.near_spec{iC}.freqs, hfr_lfr_ppc, 'black');
                % show vertical lines at points of significance
                if ~this_sig{2} == 0 
                    p22 = xline(ax5, this_sig{2}, 'LineStyle', '--', 'Color',c2);
                    p22.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{3} == 0
                    p23 = xline(ax5, this_sig{3}, 'LineStyle', '--', 'Color',c1);
                    p23.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                ax5.Box = 'off';
                ax5.YTick = [];
                ax5.XLim = [x1 100];
                ax5.XTick = [x1 10 25 50 75 100];
                ax5.Title.String = 'PPC diff';
                ax5.XLabel.String = 'Frequency (Hz)';
                ax5.Title.FontSize = 25;
                ax5.XAxis.FontSize = 18;
                ax5.XAxis.FontWeight = 'normal';
                ax5.TickDir = 'out';
                leg5 = legend({'Mean +/- 3*SD','HFR - LFR'}, 'Location','best');

                o_name = fsi_labels{iC};
                if this_sig{2} + this_sig{3} > 0
                    o_name = strcat('SIG_FSI_', o_name);
                else
                    o_name = strcat('FSI_', o_name);
                end
                WriteFig(fig, o_name, 1);
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
                this_sig = strcmp(sig_msn(:,1),msn_labels{iC});
                this_sig = sig_msn(this_sig,:);
                x1 = find(round(od.msn_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                
                % Plot STA
                fig = figure('WindowState', 'maximized');
                ax1 = subplot(2,3,1);
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
                ax2 = subplot(2,3,2);
                hold on;
                p4 = plot(ax2, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_lfr_spec{iC}.sts_vals, 'Color', c1);
                p5 = plot(ax2, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_hfr_spec{iC}.sts_vals, 'Color', c2);
                p6 = plot(ax2, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_spec{iC}.sts_vals, 'Color', c3);
                % show vertical lines at points of significance
                if ~this_sig{4} == 0 
                    p7 = xline(ax2, this_sig{4}, 'LineStyle', '--', 'Color',c2);
                    p7.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{5} == 0
                    p8 = xline(ax2, this_sig{5}, 'LineStyle', '--', 'Color',c1);
                    p8.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
                leg2 = legend({'LFR','HFR','All Trials'}, 'Location','best');
                
                % Plot PPC
                ax3 = subplot(2,3,3);
                hold on;
                p9 = plot(ax3, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_lfr_spec{iC}.ppc', 'Color', c1);
                p10 = plot(ax3, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_hfr_spec{iC}.ppc', 'Color', c2);
                p11 = plot(ax3, od.msn_res.near_spec{iC}.freqs, ...
                    od.msn_res.near_spec{iC}.ppc', 'Color', c3);
                % show vertical lines at points of significance
                if ~this_sig{2} == 0 
                    p12 = xline(ax3, this_sig{2}, 'LineStyle', '--', 'Color',c2);
                    p12.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{3} == 0
                    p13 = xline(ax3, this_sig{3}, 'LineStyle', '--', 'Color',c1);
                    p13.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
                leg3 = legend({'LFR','HFR','All Trials'}, 'Location','best');

                % Plot STS dif
                ax4 = subplot(2,3,5);
                hold on;
                p1_p2_sts = od.msn_res.near_p1_spec{iC}.sts - od.msn_res.near_p2_spec{iC}.sts;
                mean_sts = mean(p1_p2_sts,1);
                sd_sts = std(p1_p2_sts);
                hfr_lfr_sts = od.msn_res.near_hfr_spec{iC}.sts_vals - od.msn_res.near_lfr_spec{iC}.sts_vals;
                p14 = plot(ax4, od.msn_res.near_spec{iC}.freqs, mean_sts + (3*sd_sts), '--m');
                p15 = plot(ax4, od.msn_res.near_spec{iC}.freqs, mean_sts - (3*sd_sts), '--m');
                p15.Annotation.LegendInformation.IconDisplayStyle = 'off';
                p16 = plot(ax4, od.msn_res.near_spec{iC}.freqs, hfr_lfr_sts, 'black');
                % show vertical lines at points of significance
                if ~this_sig{4} == 0 
                    p17 = xline(ax4, this_sig{4}, 'LineStyle', '--', 'Color',c2);
                    p17.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{5} == 0
                    p18 = xline(ax4, this_sig{5}, 'LineStyle', '--', 'Color',c1);
                    p18.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end 
                ax4.Box = 'off';
                ax4.YTick = [];
                ax4.XLim = [x1 100];
                ax4.XTick = [x1 10 25 50 75 100];
                ax4.Title.String = 'STS diff';
                ax4.XLabel.String = 'Frequency (Hz)';
                ax4.Title.FontSize = 25;
                ax4.XAxis.FontSize = 18;
                ax4.XAxis.FontWeight = 'normal';
                ax4.TickDir = 'out';
                leg4 = legend({'Mean +/- 3*SD','HFR - LFR'}, 'Location','best');


                % Plot PPC dif
                ax5 = subplot(2,3,6);
                hold on;
                p1_p2_ppc = od.msn_res.near_p1_spec{iC}.ppc - od.msn_res.near_p2_spec{iC}.ppc;
                mean_ppc = mean(p1_p2_ppc, 1);
                sd_ppc = std(p1_p2_ppc);
                hfr_lfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc' - od.msn_res.near_lfr_spec{iC}.ppc';
                p19 = plot(ax5, od.msn_res.near_spec{iC}.freqs, mean_ppc + (3*sd_ppc), '--m');
                p20 = plot(ax5, od.msn_res.near_spec{iC}.freqs, mean_ppc - (3*sd_ppc), '--m');
                p20.Annotation.LegendInformation.IconDisplayStyle = 'off';
                p21 = plot(ax5, od.msn_res.near_spec{iC}.freqs, hfr_lfr_ppc, 'black');
                % show vertical lines at points of significance
                if ~this_sig{2} == 0 
                    p22 = xline(ax5, this_sig{2}, 'LineStyle', '--', 'Color',c2);
                    p22.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                if ~this_sig{3} == 0
                    p23 = xline(ax5, this_sig{3}, 'LineStyle', '--', 'Color',c1);
                    p23.Annotation.LegendInformation.IconDisplayStyle = 'off';
                end
                ax5.Box = 'off';
                ax5.YTick = [];
                ax5.XLim = [x1 100];
                ax5.XTick = [x1 10 25 50 75 100];
                ax5.Title.String = 'PPC diff';
                ax5.XLabel.String = 'Frequency (Hz)';
                ax5.Title.FontSize = 25;
                ax5.XAxis.FontSize = 18;
                ax5.XAxis.FontWeight = 'normal';
                ax5.TickDir = 'out';
                leg5 = legend({'Mean +/- 3*SD','HFR - LFR'}, 'Location','best');

                o_name = msn_labels{iC};
                if this_sig{2} + this_sig{3} > 0
                    o_name = strcat('SIG_MSN_', o_name);
                else
                    o_name = strcat('MSN_', o_name);
                end
                WriteFig(fig, o_name, 1);
                close;
            end
        end
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
