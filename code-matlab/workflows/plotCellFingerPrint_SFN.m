cd('E:\Dropbox (Dartmouth College)\AnalysisResults\FieldTripResults\ft_results');

fw = [5,100];

set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesFontWeight','bold')



sig_msn_label = 'R119-2007-07-10-TT11_3';
sig_fsi_label = 'R117-2007-06-18-TT07_2';
nonsig_msn_label = 'R132-2007-10-19-TT01_4';
nonsig_fsi_label = 'R117-2007-06-01-TT07_1';

%% Uncomment and run only if labels are changed
% iC = 0;
% sfn_msn_data = {[],[]}; %1 for sig, 2 for nonsig
% sfn_fsi_data = {[],[]}; %1 for sig, 2 for nonsig
% rats = {'R117','R119','R131','R132'};
% pk_thresh = -1;
% load('cellsOfInterest_cw.mat');
% load('significant_dif.mat');

% for idx = 1:length(rats)
%     curRat = rats{idx};
%     searchString = strcat(curRat,'*ft_spec.mat');
%     ofiles = dir(searchString);
%     for jdx = 1:length(ofiles)        
%         load(ofiles(jdx).name);
%         fsi_labels  = od.label(od.cell_type == 2);
%         for iF = 1:length(fsi_labels)
%             this_label = extractBefore(fsi_labels{iF}, '.t');
%             if strcmp(this_label, sig_fsi_label)
%                 iC = 1;
%             elseif strcmp(this_label, nonsig_fsi_label)
%                 iC = 2;
%             else
%                 continue;
%             end
%             
%             sfn_fsi_data{iC}.hfr_sta = od.fsi_res.near_hfr_spec{iF}.sta_vals;
%             sfn_fsi_data{iC}.lfr_sta = od.fsi_res.near_lfr_spec{iF}.sta_vals;
%             sfn_fsi_data{iC}.sta_time =  od.fsi_res.near_spec{iF}.sta_time;
%             sfn_fsi_data{iC}.hfr_spk_count = od.fsi_res.near_hfr_spec{iF}.spk_count;
%             sfn_fsi_data{iC}.lfr_spk_count = od.fsi_res.near_lfr_spec{iF}.spk_count;
%             
%             % Find PPC_peaks  
%             lf = find(od.fsi_res.near_spec{iF}.freqs >= fw(1), 1, 'first');
%             rf = find(od.fsi_res.near_spec{iF}.freqs <= fw(2), 1, 'last');
%             
%             % All trials Peak
%             pks = findpeaks(od.fsi_res.near_spec{iF}.subsampled_ppc(lf:rf)); 
%             [~, loc1] = max(od.fsi_res.near_spec{iF}.subsampled_ppc(lf+pks.loc-1));
%             near_ppc_pk = od.fsi_res.near_spec{iF}.freqs(pks.loc(loc1) + fw(1) - 1);
%             
%             % HFR Peak
%             pks = findpeaks(od.fsi_res.near_hfr_spec{iF}.subsampled_ppc(lf:rf)); 
%             [~, loc1] = max(od.fsi_res.near_hfr_spec{iF}.subsampled_ppc(lf+pks.loc-1));
%             near_hfr_ppc_pk = od.fsi_res.near_spec{iF}.freqs(pks.loc(loc1) + fw(1) - 1);
% 
%             % LFR Peak
%             pks = findpeaks(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc(lf:rf));
%             max_val = max(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc);
%             [~, loc1] = max(od.fsi_res.near_lfr_spec{iF}.subsampled_ppc(lf+pks.loc-1));
%             near_lfr_ppc_pk = od.fsi_res.near_spec{iF}.freqs(pks.loc(loc1) + fw(1) - 1);
%             
%             
%             sfn_fsi_data{iC}.hfr_ppc = od.fsi_res.near_hfr_spec{iF}.subsampled_ppc;
%             sfn_fsi_data{iC}.lfr_ppc = od.fsi_res.near_lfr_spec{iF}.subsampled_ppc;
%             sfn_fsi_data{iC}.all_ppc = od.fsi_res.near_spec{iF}.subsampled_ppc;
%             sfn_fsi_data{iC}.freqs = od.fsi_res.near_spec{iF}.freqs;
%             sfn_fsi_data{iC}.hfr_ppc_peak = near_hfr_ppc_pk;
%             sfn_fsi_data{iC}.lfr_ppc_peak = near_lfr_ppc_pk;
%             sfn_fsi_data{iC}.all_ppc_peak = near_ppc_pk;
%             
%             p1_p2_ppc = od.fsi_res.near_p1_spec{iF}.subsampled_ppc - od.fsi_res.near_p2_spec{iF}.subsampled_ppc;
%             mean_ppc = mean(p1_p2_ppc, 1);
%             sd_ppc = std(p1_p2_ppc, 1);  
%             sfn_fsi_data{iC}.mean_ppc = mean_ppc;
%             sfn_fsi_data{iC}.sd_ppc = sd_ppc;
%             sfn_fsi_data{iC}.diff_ppc = od.fsi_res.near_hfr_spec{iF}.subsampled_ppc - od.fsi_res.near_lfr_spec{iF}.subsampled_ppc;
%             
%             % Hacky and hardcoded
%             if iC == 1
%                  sfn_fsi_data{iC}.sig_freq = sfn_fsi_data{iC}.freqs(33);
%             else
%                  sfn_fsi_data{iC}.sig_freq = nan;
%             end
%         end
%         
%                 msn_labels  = od.label(od.cell_type == 1);
%         for iM = 1:length(msn_labels)
%             this_label = extractBefore(msn_labels{iM}, '.t');
%             if strcmp(this_label, sig_msn_label)
%                 iC = 1;
%             elseif strcmp(this_label, nonsig_msn_label)
%                 iC = 2;
%             else
%                 continue;
%             end
%             
%             sfn_msn_data{iC}.hfr_sta = od.msn_res.near_hfr_spec{iM}.sta_vals;
%             sfn_msn_data{iC}.lfr_sta = od.msn_res.near_lfr_spec{iM}.sta_vals;
%             sfn_msn_data{iC}.sta_time =  od.msn_res.near_spec{iM}.sta_time;
%             sfn_msn_data{iC}.hfr_spk_count = od.msn_res.near_hfr_spec{iM}.spk_count;
%             sfn_msn_data{iC}.lfr_spk_count = od.msn_res.near_lfr_spec{iM}.spk_count;
%             
%             % Find PPC_peaks  
%             lf = find(od.msn_res.near_spec{iM}.freqs >= fw(1), 1, 'first');
%             rf = find(od.msn_res.near_spec{iM}.freqs <= fw(2), 1, 'last');
%             
%             % All trials Peak
%             pks = findpeaks(od.msn_res.near_spec{iM}.ppc(lf:rf)); 
%             [~, loc1] = max(od.msn_res.near_spec{iM}.ppc(lf+pks.loc-1));
%             near_ppc_pk = od.msn_res.near_spec{iM}.freqs(pks.loc(loc1) + fw(1) - 1);
%             
%             % HFR Peak
%             pks = findpeaks(od.msn_res.near_hfr_spec{iM}.ppc(lf:rf)); 
%             [~, loc1] = max(od.msn_res.near_hfr_spec{iM}.ppc(lf+pks.loc-1));
%             near_hfr_ppc_pk = od.msn_res.near_spec{iM}.freqs(pks.loc(loc1) + fw(1) - 1);
% 
%             % LFR Peak
%             pks = findpeaks(od.msn_res.near_lfr_spec{iM}.ppc(lf:rf));
%             max_val = max(od.msn_res.near_lfr_spec{iM}.ppc);
%             [~, loc1] = max(od.msn_res.near_lfr_spec{iM}.ppc(lf+pks.loc-1));
%             near_lfr_ppc_pk = od.msn_res.near_spec{iM}.freqs(pks.loc(loc1) + fw(1) - 1);
%             
%             
%             sfn_msn_data{iC}.hfr_ppc = od.msn_res.near_hfr_spec{iM}.ppc;
%             sfn_msn_data{iC}.lfr_ppc = od.msn_res.near_lfr_spec{iM}.ppc;
%             sfn_msn_data{iC}.all_ppc = od.msn_res.near_spec{iM}.ppc;
%             sfn_msn_data{iC}.freqs = od.msn_res.near_spec{iM}.freqs;
%             sfn_msn_data{iC}.hfr_ppc_peak = near_hfr_ppc_pk;
%             sfn_msn_data{iC}.lfr_ppc_peak = near_lfr_ppc_pk;
%             sfn_msn_data{iC}.all_ppc_peak = near_ppc_pk;
%             
%             p1_p2_ppc = od.msn_res.near_p1_spec{iM}.ppc - od.msn_res.near_p2_spec{iM}.ppc;
%             mean_ppc = mean(p1_p2_ppc, 1);
%             sd_ppc = std(p1_p2_ppc, 1);  
%             sfn_msn_data{iC}.mean_ppc = mean_ppc;
%             sfn_msn_data{iC}.sd_ppc = sd_ppc;
%             sfn_msn_data{iC}.diff_ppc = od.msn_res.near_hfr_spec{iM}.ppc - od.msn_res.near_lfr_spec{iM}.ppc;
%             
%             % Hacky and hardcoded
%             if iC == 1
%                  sfn_msn_data{iC}.sig_freq = sfn_msn_data{iC}.freqs(64);
%             else
%                  sfn_msn_data{iC}.sig_freq = nan;
%             end
%         end
% 
%     end
% end
%% Do plotting stuff
load('./sfn_data.mat')

c2 = [26/255 255/255 26/255]; % Green
c1 = [75/255 0/255 146/255];  % Violet/Purple
c3 = [0.7 0.7 0.7]; %Gray
%% Plot all STA-PPC
close all;
figure
fs1 = 28;
fs2 = 28;
ax = subplot(6,2,[1 3 5]);
hold on
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_msn_data{1}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{1}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{1}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100];
ax.XTickLabel{2} = '';
ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'northeast';

ax = subplot(6,2,[9 11]);
hold on
plot(sfn_msn_data{1}.sta_time, sfn_msn_data{1}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{1}.sta_time, sfn_msn_data{1}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5];
ax.XTickLabel{2} = '';
ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';

ax = subplot(6,2,[2 4 6]);
hold on
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_msn_data{2}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{2}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{2}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'northeast';

ax = subplot(6,2,[10 12]);
hold on
plot(sfn_msn_data{2}.sta_time, sfn_msn_data{2}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{2}.sta_time, sfn_msn_data{2}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';

fig = figure(2);

ax = subplot(6,2,[1 3 5]);
hold on
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_fsi_data{1}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{1}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{1}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'north';

ax = subplot(6,2,[9 11]);
hold on
plot(sfn_fsi_data{1}.sta_time, sfn_fsi_data{1}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{1}.sta_time, sfn_fsi_data{1}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';

ax = subplot(6,2,[2 4 6]);
hold on
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_fsi_data{2}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{2}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{2}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'north';

ax = subplot(6,2,[10 12]);
hold on
plot(sfn_fsi_data{2}.sta_time, sfn_fsi_data{2}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{2}.sta_time, sfn_fsi_data{2}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';



%% Plot MSN STA-PPC
close all;
figure
fs1 = 20;
fs2 = 20;
ax = subplot(3,2,[1 3]);
hold on
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_msn_data{1}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{1}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{1}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'north';

ax = subplot(3,2,5);
hold on
plot(sfn_msn_data{1}.sta_time, sfn_msn_data{1}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{1}.sta_time, sfn_msn_data{1}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';

ax = subplot(3,2,[2 4]);
hold on
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_msn_data{2}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{2}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{2}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'north';

ax = subplot(3,2,6);
hold on
plot(sfn_msn_data{2}.sta_time, sfn_msn_data{2}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_msn_data{2}.sta_time, sfn_msn_data{2}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';

%% Plot FSI STA-PPC
close all;
figure
fs1 = 20;
fs2 = 20;
ax = subplot(3,2,[1 3]);
hold on
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_fsi_data{1}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{1}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{1}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'north';

ax = subplot(3,2,5);
hold on
plot(sfn_fsi_data{1}.sta_time, sfn_fsi_data{1}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{1}.sta_time, sfn_fsi_data{1}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';

ax = subplot(3,2,[2 4]);
hold on
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.lfr_ppc, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.hfr_ppc, 'Color', c2, 'LineWidth', 3);
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.all_ppc, 'Color', c3, 'LineWidth', 3);
q0 = xline(sfn_fsi_data{2}.lfr_ppc_peak, 'Color', c1, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{2}.hfr_ppc_peak, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{2}.all_ppc_peak, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.005 0.035];
ax.YTick = [0 0.035];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.YAxis.FontWeight = 'normal';
ax.XAxis.FontWeight = 'normal';
ax.TickDir = 'out';
ax.YAxis.Exponent = 0;
leg = legend({'LFR trials', 'HFR trials', 'ALL trials'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'north';

ax = subplot(3,2,6);
hold on
plot(sfn_fsi_data{2}.sta_time, sfn_fsi_data{2}.lfr_sta, 'Color', c1, 'LineWidth', 3);
plot(sfn_fsi_data{2}.sta_time, sfn_fsi_data{2}.hfr_sta, 'Color', c2, 'LineWidth', 3);
ax.Box = 'off';
ax.XTick = [-0.5 -0.25 0 0.25 0.5]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.YTick = {};
ax.YLabel.String = 'STA';
ax.XLabel.String = 'Time (sec)';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';

%% Plot PPC diff
close all;
figure
fs1 = 20;
fs2 = 20;

% Plot Sig MSN
ax = subplot(2,2,1);
hold on;
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.diff_ppc, 'Color', 'magenta');
plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.mean_ppc + 3*sfn_msn_data{1}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0 = plot(sfn_msn_data{1}.freqs, sfn_msn_data{1}.mean_ppc - 3*sfn_msn_data{1}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_msn_data{1}.sig_freq, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.06 0.06];
ax.YTick = [-0.06 0 0.06];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC Difference';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';
leg = legend({'HFR - LFR', 'Mean +/- 3*SD'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'south';

% Plot Nongsig MSN
ax = subplot(2,2,2);
hold on;
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.diff_ppc, 'Color', 'magenta');
plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.mean_ppc + 3*sfn_msn_data{2}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0 = plot(sfn_msn_data{2}.freqs, sfn_msn_data{2}.mean_ppc - 3*sfn_msn_data{2}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.06 0.06];
ax.YTick = [-0.06 0 0.06];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC Difference';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';
leg = legend({'HFR - LFR', 'Mean +/- 3*SD'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'south';

% Plot Sig FSI
ax = subplot(2,2,3);
hold on;
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.diff_ppc, 'Color', 'magenta');
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.mean_ppc + 3*sfn_fsi_data{1}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0 = plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.mean_ppc - 3*sfn_fsi_data{1}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
q0 = xline(sfn_fsi_data{1}.sig_freq, 'Color', c2, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.01 0.01];
ax.YTick = [-0.01 0 0.01];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC Difference';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';
leg = legend({'HFR - LFR', 'Mean +/- 3*SD'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'south';

% Plot Nongsig FSI
ax = subplot(2,2,4);
hold on;
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.diff_ppc, 'Color', 'magenta');
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.mean_ppc + 3*sfn_fsi_data{2}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0 = plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.mean_ppc - 3*sfn_fsi_data{2}.sd_ppc, 'Color', c3, 'LineStyle','--','LineWidth', 3);
q0.Annotation.LegendInformation.IconDisplayStyle = 'off';
ax.Box = 'off';
ax.XLim = fw;
ax.YAxis.Exponent = 0;
ax.YLim = [-0.01 0.01];
ax.YTick = [-0.01 0 0.01];
ax.XTick = [5 25 50 75 100]; ax.XTickLabel{2} = ''; ax.XTickLabel{4} = '';
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'PPC Difference';
ax.XAxis.FontSize = fs1;
ax.YAxis.FontSize = fs1;
ax.XAxis.FontWeight = 'normal';
ax.YAxis.FontWeight = 'normal';
ax.TickDir = 'out';
leg = legend({'HFR - LFR', 'Mean +/- 3*SD'});
leg.Box = 'off';
leg.FontSize = fs2;
leg.FontWeight = 'normal';
leg.Location = 'south';

%%


ax3.XLabel.String = 'Frequency (in Hz)';
            ax3.Title.FontSize = 25;
            ax3.XAxis.FontSize = 18;
            ax3.XLim = fw;
            ax3.XTick = [5 10 20 30 40 50 60 70 80 90 100];
            ax3.XAxis.FontWeight = 'normal';
            ax3.YAxis.FontWeight = 'normal';
            ax3.TickDir = 'out';



dummy =1 ;

subplot(2,2,2)
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.diff_ppc);
subplot(2,2,3)
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.diff_ppc);
subplot(2,2,4)
plot(sfn_fsi_data{1}.freqs, sfn_fsi_data{1}.diff_ppc);


%%
dummy = 1;
% pbaspect([1 0.3 1])
subplot(3,2,5)
plot(sfn_fsi_data{1}.sta_time, sfn_fsi_data{1}.hfr_sta);
% pbaspect([1 0.3 1])
subplot(3,2,[2 4])
plot(sfn_fsi_data{2}.freqs, sfn_fsi_data{2}.hfr_ppc);
% pbaspect([1 0.3 1])
subplot(3,2,6)
plot(sfn_fsi_data{2}.sta_time, sfn_fsi_data{2}.hfr_sta);
% pbaspect([1 0.3 1])

