function AMPX_plane_count_bar(plane_stats, plane_stats_control, band_name)
%% plot some stuff
F = figure(200);
set(gcf, 'PaperPositionMode', 'auto', 'color', 'w')
set(F, 'Position', [200, 200, 900 700])

% y = [ mean(low_gamma.control) mean(high_gamma.control); mean(low_gamma.ipsi) mean(high_gamma.ipsi); mean(low_gamma.contra) mean(high_gamma.contra)];
% E = [low_control_SEM, high_control_SEM ; low_ipsi_SEM, high_ipsi_SEM; low_contra_SEM, high_contra_SEM];
y = [mean(plane_stats.rsq2); mean(plane_stats_control.rsq2)];
E = [std(plane_stats.rsq2); std(plane_stats_control.rsq2)];
colors = linspecer(4);
b= boxplot([plane_stats.rsq2; plane_stats_control.rsq2]', 'labels', {band_name, 'Random'}, 'notch', 'on');
% b = b(y, .6);
h=findobj(gca,'tag','Outliers');
set(h,'Visible','off')

% set(b(6,1), 'color', 'k', 'linewidth', 2)
% set(b(6,2), 'color', 'k', 'linewidth', 2)
set(b(1:6,1), 'color', 'k', 'linewidth', 2)
set(b(1:6,2), 'color', colors(2,:), 'linewidth', 2)
set(findobj(get(h(1), 'parent'), 'type', 'text'), 'fontname', 'helvetica','fontsize', 16)% 'fontweight', 'demi', 'position', [1 1 1 1]);
% set(b(2), 'facecolor', colors(2,:))
xlabel('Event Type','FontSize', 16, 'fontweight', 'demi','fontname', 'helvetica')
ylim([0 100])
% set(b(1), 'facecolor', colors(2,:))
% set(b(2), 'facecolor', colors(3,:))
% hold on
% h = errorbar(y, E);
% for i  = 1
% set(h(i),'color','k', 'LineStyle','none', 'linewidth', 2)
% end
d_x = get(h(1), 'xData');
% d_x2 = get(h(2), 'xData');
% set(h(1), 'xData', d_x-0.15)
% set(h(2), 'xData', d_x+0.15)
ylabel({'Percentage of varience explained  '; 'by r^2 Across for Gamma Events' }, 'fontsize', 16)
set(gca, 'fontsize', 16, 'FontName', 'helvetica','ytick', [0:10:100])
set(gcf, 'color', 'w')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
%%
a= plane_stats.rsq2;
b =plane_stats_control.rsq2;
[h, p, ci, stats] = ttest2(plane_stats.rsq2, plane_stats_control.rsq2);
fprintf(['\n' band_name ' Gamma\nmean=' num2str(nanmean(a)) ', SD±' num2str(std(a)) '; independent samples t-test: t(' num2str(stats.df) ')' ' = ' num2str(stats.tstat), ', p= ' num2str(p) '\n'])

fprintf(['\nRandom \nmean=' num2str(nanmean(b)) ', SD±' num2str(std(b)) '\n'])
figure(200)
hold on
l_width = 2;
% if h == 1;
%     ip_con = 1.7:0.001:1.85;
%     ip_con(2,:) = ones(1,length(ip_con));
%     plot(ip_con(2,:).*2, ip_con(1,:), '-k','linewidth', l_width);
%     plot((ip_con(2,:).*3), ip_con(1,:), '-k','linewidth', l_width)
%     plot(2:3, [1.85 1.85], '-k', 'linewidth', l_width)
% end


tx_x = 3.3;
text(tx_x, -.155, '** {\itp}<0.005', 'FontSize', 14)
text(tx_x, -.1, ' *  {\itp}<0.05', 'FontSize', 14)
text(tx_x-0.15, -.2, 'two sample t-test', 'FontSize', 14)


%%
saveas(F, ['G:\Naris\paper_figs\SI_plane_count.fig'])
print(F, '-dpng','-r300',['G:\Naris\paper_figs\SI_plane_count.png'])