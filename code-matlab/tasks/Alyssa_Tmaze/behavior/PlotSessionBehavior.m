function PlotSessionBehavior(cfg_in,metadata,ExpKeys)
% function PlotSessionBehavior(cfg_in,metadata,ExpKeys)
%

cfg_def.ax = [];
cfg_def.order = {'pedL','pedR','trialL','trialR'};
cfg_def.cols = {[0.7 0.1 0.1],[0.1 0.1 0.7],[1 0 0],[0 0 1]};
cfg_def.lw = 2;
cfg_def.fs = 14;

cfg = ProcessConfig2(cfg_def,cfg_in);

%
pedL_idx = find(strcmp('pedL',cfg.order));
pedR_idx = find(strcmp('pedR',cfg.order));
trialL_idx = find(strcmp('trialL',cfg.order));
trialR_idx = find(strcmp('trialR',cfg.order));

%
if ~isempty(cfg.ax)
    axes(cfg.ax);
else
   figure; 
   subplot(311);
end

% pre/post
switch ExpKeys.Pedestal
    case 'L'
        this_idx = pedL_idx;
    case 'R'
        this_idx = pedR_idx;
end

plot([ExpKeys.prerecord(1) ExpKeys.prerecord(2)],[this_idx this_idx],'LineWidth',cfg.lw,'Color',cfg.cols{this_idx});
hold on;
plot([ExpKeys.postrecord(1) ExpKeys.postrecord(2)],[this_idx this_idx],'LineWidth',cfg.lw,'Color',cfg.cols{this_idx});

% trials
for iT = 1:length(metadata.taskvars.sequence)
    
    switch metadata.taskvars.sequence{iT}
        case 'L'
            this_idx = trialL_idx;
        case 'R'
            this_idx = trialR_idx;
    end
    
    plot([metadata.taskvars.trial_iv.tstart(iT) metadata.taskvars.trial_iv.tend(iT)],[this_idx this_idx],'LineWidth',cfg.lw,'Color',cfg.cols{this_idx});
    
end

% pedestals
for iT = 1:length(metadata.taskvars.rest_iv_pedL.tstart)
    
    this_idx = pedL_idx;
    
    plot([metadata.taskvars.rest_iv_pedL.tstart(iT) metadata.taskvars.rest_iv_pedL.tend(iT)],[this_idx this_idx],'LineWidth',cfg.lw,'Color',cfg.cols{this_idx});
    
end

for iT = 1:length(metadata.taskvars.rest_iv_pedR.tstart)
    
    this_idx = pedR_idx;
    
    plot([metadata.taskvars.rest_iv_pedR.tstart(iT) metadata.taskvars.rest_iv_pedR.tend(iT)],[this_idx this_idx],'LineWidth',cfg.lw,'Color',cfg.cols{this_idx});
    
end

% set
set(gca,'YLim',[0.5 length(cfg.order)+0.5],'FontSize',cfg.fs,'YTick',1:length(cfg.order),'YTickLabel',cfg.order, ...
    'XLim',[ExpKeys.prerecord(1) ExpKeys.postrecord(2)]);
box off; xlabel('time(s)');

title(sprintf('%s (%s-restricted)',ExpKeys.goodSWR{1}(1:15),ExpKeys.RestrictionType));