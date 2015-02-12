function tch = TCPlot(cfg_in,TC)
% function TCPlot(cfg,TC)
%
% plot tuning curves, ordered by peak location
%
% INPUTS
%
% TC: tuning curves object (as returned by MakeTC())
%
% OUTPUTS
%
% tch: vector of handles to plotted tuning curves (so can do things like,
%  set(tch,'FaceColor',[1 0 0])
%
% CONFIGS
%
% cfg_def.ax = []; % if specified, plot here; otherwise create new figure
%
% disclaimer:
% 
% MvdM rough version copied from youkitan workflow chunk; many improvements
% to make especially options to plot all tuning curves, only those with
% fields, etc..

cfg_def = [];
cfg_def.ax = []; % if specified, plot here; otherwise create new figure
cfg = ProcessConfig2(cfg_def,cfg_in);

if isempty(cfg.ax)
    figure;
else
    axes(cfg.ax);
end

tc_temp = TC.tc(:,TC.template_idx)';
clear tch;
for iC = size(tc_temp,1):-1:1       
    subaxis(size(tc_temp,1),1,iC,'SpacingVert',0);
    tch(iC) = area(tc_temp(iC,:)); hold on; set(tch(iC),'FaceColor',[0 0 0]);
    ylims = get(gca,'Ylim');
    plot([TC.peak_idx(iC) TC.peak_idx(iC)],[ylims(1) ylims(2)],'r:');
    set(gca, 'XTick', [],'YTick',[]); 
    yl = ylabel(iC); 
    set(yl,'Rotation',0,'Fontsize',8);
    %if iC == 1; title(ENC_data(iT).trial_type); end
end
x = size(tc_temp,2);
%set(gca,'XTick',[1 ceil(x/2) x],'XTicklabel',[1 137 274],'Ticklength', [0 0],'Fontsize',8);


