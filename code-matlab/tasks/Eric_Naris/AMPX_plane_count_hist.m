function [h2, stats] = AMPX_plane_count_hist(plane_stats, plane_stats_control, band_name)
a = plane_stats.rsq2;
b = plane_stats_control.rsq2;

[h, p, ci, stats] = ttest2(a, b)
figure
[n1, p1] = hist(b, 50);
[n2, n] = hist(a, 50);
h = bar(n, n1, 'r', 'barwidth', .9, 'edgecolor', 'none');
hold on
h2 = bar(n, n2, 'b', 'barwidth', .9, 'edgecolor', 'none');
legend([h h2], {'Control', 'Gamma'})

h2 = findobj(gca,'Type','patch');
set(h2,'facealpha',0.25)
set(gca, 'fontsize', 16,'FontName', 'helvetica','xtick', [0:10:90])
xlabel({['Percentage of 2D ' band_name  ' gamma power varience']; 'explained by best fit plane'}, 'Fontsize', 16, 'fontname', 'helvetica')
% hold on
% bar(hist(b,100))
% ha = hist(a,100)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b','EdgeColor','w','facealpha',0.25)
% hold on;
% hb = hist(b,100)
% h1 = findobj(hb,'Type','patch');
% set(h1,'Facecolor', 'b', 'facealpha',0.25);

fprintf([band_name ' gamma: mean = ' num2str(mean(a)) '\n']);
fprintf(['Random ' band_name ' gamma: mean = ' num2str(mean(b)) '\n'])

fprintf(['\nLow Gamma\nmean=' num2str(nanmean(a)) ', SD±' num2str(std(a)) '; independent samples t-test: t(' num2str(stats.df) ')' ' = ' num2str(stats.tstat), ', p= ' num2str(p) '\n'])

fprintf(['\nRandom \nmean=' num2str(nanmean(b)) ', SD±' num2str(std(b)) '\n'])