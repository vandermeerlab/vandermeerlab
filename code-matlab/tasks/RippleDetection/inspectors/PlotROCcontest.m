function PlotROCcontest(cfg,ROCdata1,ROCdata2,varargin)
%PLOTROCCONTEST Summary of this function goes here
%   Detailed explanation goes here

mfun = mfilename;
cfg_def.rating = 1;
cfg = ProcessConfig(cfg_def,cfg,mfun);
col = linspecer(6);

figure; hold on; title('ROC method compare')
plot(ROCdata1.roc6,ROCdata1.(['roc',num2str(cfg.rating)]),'Color',col(1,:),'linewidth',2); name = inputname(2); if isempty(name); name = ['Input ',num2str(1)]; end; lgnd{1} = name;
plot(ROCdata2.roc6,ROCdata2.(['roc',num2str(cfg.rating)]),'Color',col(2,:),'linewidth',2); name = inputname(3); if isempty(name); name = ['Input ',num2str(2)]; end; lgnd{2} = name;
% Note that inputname(1) is cfg, so we need to use inputname(iArg+1) to get the name of the input we want

for iVarg = 1:length(varargin)
    iArg = iVarg+2;
    ROCdata = varargin{iVarg};
    name = inputname(iArg+1); if isempty(name); name = ['Input ',num2str(iArg)]; end; lgnd{iArg} = name;
    plot(ROCdata.roc6,ROCdata.(['roc',num2str(cfg.rating)]),'Color',col(3,:),'linewidth',2)
end

xlabel('False positive rate')
ylabel('True positive rate')
legend(lgnd,'location','southeast');

set(gca,'YLim',[0 1.05],'XLim',[-0.05 0.6])

end

