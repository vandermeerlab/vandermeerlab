function AMPX_Naris_fig_5_gamma_count()
%% AMPX_Naris_fig_5_gamma_count: creates a bar plot for the number of
% gamma event occurances across the "control", ipsilateral and
% contralateral phases.  requires that the count has been created by
% Naris_gamma_count
%
%
% EC - 2016-06-05


%% collect the output stats
load('G:\Naris\Paper_naris_gamma.mat')
%%
phases = {'control', 'pre', 'right', 'left', 'post'};
loop = 1;
sess_name = fieldnames(naris_gamma);
all_count_low = [];
all_count_high = [];
for iSess = 1:length(sess_name)
    
    all_count_low = [all_count_low; naris_gamma.(sess_name{iSess}).count_L];
    all_count_high = [all_count_high; naris_gamma.(sess_name{iSess}).count_H];
    
end

%% Normalize to Control
for iphase = 1:5
    n_all_count_low(:,iphase) = all_count_low(:,iphase)./all_count_low(:,1);
    n_all_count_high(:,iphase) = all_count_high(:,iphase)./all_count_high(:,1);
    low_SEM(iphase) = std(n_all_count_low(:,iphase))/(sqrt(length(all_count_low(:,iphase))));
    high_SEM(iphase) = std(n_all_count_high(:,iphase))/(sqrt(length(all_count_high(:,iphase))));
    low_mean(iphase) = mean(n_all_count_low(:,iphase));
    high_mean(iphase) = mean(n_all_count_high(:,iphase));
end
%% make some plots
F = figure(200);
set(gcf, 'PaperPositionMode', 'auto', 'color', 'w');
set(F, 'Position', [200, -50, 900 700]);
x = [.86, 1.14; 1.86, 2.14; 2.86, 3.14];
y = [low_mean(1),low_mean(3), low_mean(4); high_mean(1),high_mean(3), high_mean(4)]';
E = [low_SEM(1),low_SEM(3), low_SEM(4); high_SEM(1),high_SEM(3), high_SEM(4)]';
colors = linspecer(3);

b = bar(y,.85);
set(b(1), 'facecolor', colors(1,:))
set(b(2), 'facecolor', colors(3,:))
hold on
h = errorbar(x,y, E);
for ii  = 1:2
    set(h(ii),'color','k', 'LineStyle','none', 'linewidth', 2);
end
ylim([0 2])
ylabel('Normalized Event Count')
set(gca, 'fontsize', 18, 'fontname', 'helvetica', 'fontweight', 'demi', 'xticklabel', {'Control', 'Ipsi', 'Contra'});
set(gca, 'ytick', [0:0.5:2]);
legend(b, {'low', 'high'}, 'orientation', 'horizontal');

%% Kolomogorov-Smirnoov (test for normality in the data)
ks = []; ksh = [];
for iPhase = [1 3 4];
    labels = {'control', 'pre', 'ipsi', 'contra', 'post'};
    h = kstest(all_count_low(:,iPhase));
    if h ~=1
        disp('***************************************************************')
        disp(['KS test FAIL for low ' labels{iPhase}])
        disp('***************************************************************')
        ks = [ks ; 1];
    end
    ks = [ks; 0];

    hh = kstest(all_count_high(:,iPhase));
    if hh ~=1
        disp('***************************************************************')
        disp(['KS test FAIL for high ' sites{iPhase}])
        disp('***************************************************************')
        ksh = [ksh ; 1];
    end
    ksh = [ksh; 0];
end

%% tests for differnences using paired t-test if passes KS test, or Wilcoxin sign-rank if non-parameteric.
if sum(ks)>=1
    [h_ip_con, p_ip_con] = signrank(all_count_low(:, 3), all_count_low(:, 4));
    [h_ip_ctr, p_ip_ctr] = signrank(all_count_low(:, 3), all_count_low(:, 1));
    [h_con_ctr, p_con_ctr] = signrank(all_count_low(:, 4), all_count_low(:, 1));
    
    [h_h_ip_con, h_p_ip_con] = signrank(all_count_high(:, 3), all_count_high(:, 4));
    [h_h_ip_ctr, h_p_ip_ctr] = signrank(all_count_high(:, 3), all_count_high(:, 1));
    [h_h_con_ctr, h_p_con_ctr] = signrank(all_count_high(:, 4), all_count_high(:, 1));
else
    disp('Using T-Test')
    [h_ip_con, p_ip_con, ~ ,l_stats_ip_con] = ttest2(all_count_low(:, 3), all_count_low(:, 4));
    [h_ip_ctr, p_ip_ctr, ~ ,l_stats_ip_ctr] = ttest2(all_count_low(:, 3), all_count_low(:, 1));
    [h_con_ctr, p_con_ctr, ~ , l_stats_con_ctr] = ttest2(all_count_low(:, 4), all_count_low(:, 1));
    
    [h_h_ip_con, h_p_ip_con,~ , h_stats_ip_con] = ttest2(all_count_high(:, 3), all_count_high(:, 4));
    [h_h_ip_ctr, h_p_ip_ctr, ~, h_stats_ip_ctr] = ttest2(all_count_high(:, 3), all_count_high(:, 1));
    [h_h_con_ctr, h_p_con_ctr, ~, h_stats_con_ctr] = ttest2(all_count_high(:, 4), all_count_high(:, 1));
end

%% display the outcome of the stats test
fileID = fopen('G:\Naris\Naris_stats.txt','w');
fprintf(fileID,datestr(date, 'yyyy-mm-dd-HH'))
fprintf(fileID, ['\nGamma Count Statistics'])
if sum(ks) >=1
    fprintf(fileID, '\nWilcoxin Sign Rank test\n')
    fprintf(fileID,['Ipsilateral   vs. Contralateral:    P:' num2str(p_ip_con, '%4.4f') '\n' ])
    fprintf(fileID,['Ipsilateral   vs. Control:          P:' num2str(p_ip_ctr, '%4.4f') '\n' ])
    fprintf(fileID, ['Contralateral vs. Control:          P:' num2str(p_con_ctr, '%4.4f') '\n' ])
else
    fprintf(fileID, '\n\nPaired T-Test\n')
    fprintf(fileID, ['Low Gamma Basic Stats\n'])
    fprintf(fileID, ['Ipsilateral:      Mean ' num2str(low_mean(3), '%4.4f') '   SD +/-' num2str(std(n_all_count_low(:,3)), '%4.4f') '   SEM +/-' num2str(low_SEM(3)) '\n'])
    fprintf(fileID, ['Contralateral:    Mean ' num2str(low_mean(4), '%4.4f') '   SD +/-' num2str(std(n_all_count_low(:,4)), '%4.4f') '   SEM +/-' num2str(low_SEM(4)) '\n'])
    fprintf(fileID, ['Control:          Mean ' num2str(low_mean(1), '%4.4f') '   SD +/-' num2str(std(n_all_count_low(:,1)), '%4.4f') '   SEM +/-' num2str(low_SEM(1)) '\n'])
    fprintf(fileID, '----------------------------------\n')
    fprintf(fileID, ['Ipsilateral   vs. Contralateral:    df(' num2str(l_stats_ip_con.df) ')   t:' num2str(l_stats_ip_con.tstat, '%4.4f') '  P:' num2str(p_ip_con, '%4.4f') '\n' ])
    fprintf(fileID, ['Ipsilateral   vs. Control:          df(' num2str(l_stats_ip_ctr.df) ')   t:' num2str(l_stats_ip_ctr.tstat, '%4.4f') '  P:' num2str(p_ip_ctr, '%4.4f') '\n' ])
    fprintf(fileID, [ 'Contralateral vs. Control:         df(' num2str(l_stats_con_ctr.df) ') t:' num2str(l_stats_con_ctr.tstat, '%4.4f') '  P:' num2str(p_con_ctr, '%4.4f') '\n' ])
    
    fprintf(fileID, ['\nHigh Gamma Basic Stats\n'])
    fprintf(fileID, ['Ipsilateral:      Mean ' num2str(high_mean(3), '%4.4f') '   SD +/-' num2str(std(n_all_count_high(:,3)), '%4.4f') '   SEM +/-' num2str(high_SEM(3)) '\n'])
    fprintf(fileID, ['Contralateral:    Mean ' num2str(high_mean(4), '%4.4f') '   SD +/-' num2str(std(n_all_count_high(:,4)), '%4.4f') '   SEM +/-' num2str(high_SEM(4)) '\n'])
    fprintf(fileID, ['Control:          Mean ' num2str(high_mean(1), '%4.4f') '   SD +/-' num2str(std(n_all_count_high(:,1)), '%4.4f') '   SEM +/-' num2str(high_SEM(1)) '\n'])
        fprintf(fileID, '----------------------------------\n')
    fprintf(fileID, ['Ipsilateral   vs. Contralateral:   df(' num2str(h_stats_ip_con.df) ')   t:' num2str(h_stats_ip_con.tstat, '%4.4f') '  P:' num2str(h_p_ip_con, '%4.4f') '\n' ])
    fprintf(fileID, ['Ipsilateral   vs. Control:         df(' num2str(h_stats_ip_ctr.df) ')   t:' num2str(h_stats_ip_ctr.tstat, '%4.4f') '  P:' num2str(h_p_ip_ctr, '%4.4f') '\n' ])
    fprintf(fileID, ['Contralateral vs. Control:         df(' num2str(h_stats_con_ctr.df) ')  t:' num2str(h_stats_con_ctr.tstat, '%4.4f') '  P:' num2str(h_p_con_ctr, '%4.4f') '\n' ])
end
fclose(fileID);
type 'G:\Naris\Naris_stats.txt'
copyfile('G:\Naris\Naris_stats.txt','C:\Users\mvdmlab\Dropbox\Naris Paper June 2016\Naris_stats.txt', 'f')
%% Add bars for stats tests
l_width = 1;
% high
if h_ip_con == 1;
    ip_con = (1.35:0.001:1.45)+0.15;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:).*x(2,2), ip_con(1,:), '-k','linewidth', l_width);
    plot((ip_con(2,:).*x(3,2)), ip_con(1,:), '-k','linewidth', l_width)
    plot(x(2,2):x(3,2), [ip_con(1,end) ip_con(1,end)], '-k', 'linewidth', l_width)
    if p_ip_con <0.005; text(median([x(2,2), x(3,2)])-.075, ip_con(1,end)+0.01, '**', 'FontSize', 46); elseif p_ip_con >0.0051 && p_ip_con <0.05; text(median([x(2,2), x(3,2)])-.05, ip_con(1,end)+0.01, '*', 'FontSize', 46) ;end
end

if h_ip_ctr == 1;
    ip_ctr = (1.2:0.001:1.3)+0.15;
    ip_ctr(2,:) = ones(1,length(ip_ctr));
    plot(ip_ctr(2,:).*x(2,2), ip_ctr(1,:), '-k','linewidth', l_width);
    plot((ip_ctr(2,:).*x(1,2)), ip_ctr(1,:), '-k','linewidth', l_width)
    plot(x(1,2):x(2,2), [ip_ctr(1,end) ip_ctr(1,end)], '-k', 'linewidth', l_width)
    if h_p_ip_ctr <0.005; text(median([x(1,2), x(2,2)])-.075, ip_ctr(1,end)+0.01, '**', 'FontSize', 46); elseif p_ip_ctr >0.0051 && p_ip_ctr <0.05; text(median([x(1,2), x(2,2)])-.05, ip_ctr(1,end)+0.01, '*', 'FontSize', 46) ;end
end

if h_h_con_ctr == 1;
    con_ctr = (1.35:0.001:1.45)+0.15;
    con_ctr(2,:) = ones(1,length(con_ctr));
    plot(con_ctr(2,:).*x(1,2), con_ctr(1,:), '-k','linewidth', l_width);
    plot((con_ctr(2,:).*x(3,2)), con_ctr(1,:), '-k','linewidth', l_width)
    plot(x(1,2):x(3,2), [con_ctr(1,end) con_ctr(1,end)], '-k', 'linewidth', l_width)
    if h_p_ip_ctr <0.005; text(median([x(1,2), x(3,2)])-.075, con_ctr(1,end)+0.01, '**', 'FontSize', 46); elseif p_ip_ctr >0.0051 && p_ip_ctr <0.05; text(median([x(1,2), x(3,2)])-.05, con_ctr(1,end)+0.01, '*', 'FontSize', 46) ;end
end
% low
if h_ip_con == 1;
    ip_con = 1.2:0.001:1.3;
    ip_con(2,:) = ones(1,length(ip_con));
    plot(ip_con(2,:).*x(2,1), ip_con(1,:), '-k','linewidth', l_width);
    plot((ip_con(2,:).*x(3,1)), ip_con(1,:), '-k','linewidth', l_width)
    plot(x(2,1):x(3,1), [ip_con(1,end) ip_con(1,end)], '-k', 'linewidth', l_width)
    if p_ip_con <0.005; text(median([x(2,1), x(3,1)])-.075, ip_con(1,end)+0.01, '**', 'FontSize', 46); elseif p_ip_con >0.0051 && p_ip_con <0.05; text(median([x(2,1), x(3,1)])-.05, ip_con(1,end)+0.01, '*', 'FontSize', 46) ;end
end

if h_ip_ctr == 1;
    ip_ctr = 1.05:0.001:1.15;
    ip_ctr(2,:) = ones(1,length(ip_ctr));
    plot(ip_ctr(2,:).*x(2,1), ip_ctr(1,:), '-k','linewidth', l_width);
    plot((ip_ctr(2,:).*x(1,1)), ip_ctr(1,:), '-k','linewidth', l_width)
    plot(x(1,1):x(2,1), [ip_ctr(1,end) ip_ctr(1,end)], '-k', 'linewidth', l_width)
    if p_ip_ctr <0.005; text(median([x(1,1), x(2,1)])-.075, ip_ctr(1,end)+0.01, '**', 'FontSize', 46); elseif p_ip_ctr >0.0051 && p_ip_ctr <0.05; text(median([x(1,1), x(2,1)])-.05, ip_ctr(1,end)+0.01, '*', 'FontSize', 46) ;end
end

if h_con_ctr == 1;
    con_ctr = 1.45:0.001:1.55;
    con_ctr(2,:) = ones(1,length(con_ctr));
    plot(con_ctr(2,:).*x(1,1), con_ctr(1,:), '-k','linewidth', l_width);
    plot((con_ctr(2,:).*x(3,1)), con_ctr(1,:), '-k','linewidth', l_width)
    plot(x(1,1):x(3,1), [con_ctr(1,end) con_ctr(1,end)], '-k', 'linewidth', l_width)
    if p_ip_ctr <0.005; text(median([x(1,1), x(3,1)])-.075, con_ctr(1,end)+0.01, '**', 'FontSize', 46); elseif p_ip_ctr >0.0051 && p_ip_ctr <0.05; text(median([x(1,1), x(3,1)])-.05, con_ctr(1,end)+0.01, '*', 'FontSize', 46) ;end
end

tx_x = 3.3;
text(tx_x, -.155, '** {\itp}<0.005', 'FontSize', 14)
text(tx_x, -.1, ' *  {\itp}<0.05', 'FontSize', 14)
if sum(ks) >=1 ;text(tx_x-0.15, -.22, 'Wilcoxon rank sum', 'FontSize', 14); else text(tx_x-0.05, -.22, 'Paired t-test', 'FontSize', 14); end

%% save the plot
saveas(gcf, 'D:\DATA\Paper_figs\Fig5_D', 'fig')
saveas(gcf, 'D:\DATA\Paper_figs\Fig5_D', 'epsc')
