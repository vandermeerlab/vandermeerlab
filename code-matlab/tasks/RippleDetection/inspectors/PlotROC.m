function PlotROC(cfg,ROCdata)
%PlotROC Plot ROC data
%   Detailed explanation goes here

% get area under the curve to generate a score
% trapz(), cumsum(), cumtrapz()...(lol)

score(1) = trapz(ROCdata.roc1);
score(2) = trapz(ROCdata.roc2);
score(3) = trapz(ROCdata.roc3);
score(4) = trapz(ROCdata.roc4);
score(5) = trapz(ROCdata.roc5);

cfg_temp = []; cfg_temp.target = 'GetROCdata'; cfg_temp.parameter = 'thresholds'; cfg_temp.verbose = 0;
thresholds = GetHistory(cfg_temp,ROCdata); threshold_diff = max(thresholds{1,1}) - min(thresholds{1,1});

ROCscore = sum(score)/threshold_diff;

col = linspecer(6);

figure; hold on; title(sprintf('ROC curve, score %.2f',ROCscore),'FontSize',14)

plot(ROCdata.roc6,ROCdata.roc1,'Color',col(1,:),'linewidth',2)
plot(ROCdata.roc6,ROCdata.roc2,'Color',col(2,:),'linewidth',2)
plot(ROCdata.roc6,ROCdata.roc3,'Color',col(3,:),'linewidth',2)
plot(ROCdata.roc6,ROCdata.roc4,'Color',col(4,:),'linewidth',2)
plot(ROCdata.roc6,ROCdata.roc5,'Color',col(5,:),'linewidth',2)

xlabel('False positive rate')
ylabel('True positive rate')
legend('1','2','3','4','5','location','southeast');

set(gca,'YLim',[0 1.05],'XLim',[-0.05 0.6])

end

