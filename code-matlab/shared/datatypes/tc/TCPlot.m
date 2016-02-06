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
% cfg_def.order = 1; % if 1, order by template_idx; 2, order by field_template_idx
% cfg_def.mode = 'area'; % 'line', 'area', 'heat'
% cfg_def.color = [0 0 0];
% cfg_def.alpha = 0;
% cfg_def.binsize = 1; % how many cm in each bin
% cfg_def.cp = []; % location of choice point
% cfg_def.YOrd = 'normal'; % order of Y axis plotting, 'normal' in 
%      decreasing order or 'reverse' in increasing order. 
% cfg_def.XDir = 'normal'; % 'normal', 'reverse' to flip horizontal
%
% disclaimer:
% 
% MvdM rough version copied from youkitan workflow chunk; many improvements
% to make especially options to plot all tuning curves, only those with
% fields, etc..
% aacarey edit sept 2015 for XDir and YOrd

cfg_def = [];
cfg_def.ax = []; % if specified, plot here; otherwise create new figure
cfg_def.order = 1; % if 1, order by template_idx; 2, order by field_template_idx
cfg_def.mode = 'area'; % 'line'
cfg_def.color = [0 0 0];
cfg_def.alpha = 0;
cfg_def.binsize = 1;
cfg_def.cp = [];
cfg_def.YOrd = 'normal'; % order of plotting along Y, 'normal' in decreasing order or 'reverse' in increasing order
cfg_def.XDir = 'normal'; % 'normal' 'reverse'

cfg = ProcessConfig2(cfg_def,cfg_in);

%
if isempty(cfg.ax)
    figure;
else
    axes(cfg.ax);
end

%
switch cfg.order
    case 1
        tc_temp = TC.tc(:,TC.template_idx)';
    case 2
        tc_temp = TC.tc(:,TC.field_template_idx)';
    otherwise
        tc_temp = TC.tc';
end

switch cfg.YOrd % handle reverse Y axis plotting of tuning curves, if specified (aacarey)
    case 'normal'
        cellNum = 1:size(tc_temp,1);
        cellNum = strtrim(cellstr(num2str(cellNum'))'); % convert to cell array of strings cellNum = {'3' '2' '1'}
    case 'reverse'
        cellNum = size(tc_temp,1):-1:1; % we flip the cell numbering and plot order for cell 1 to be at bottom
        cellNum = strtrim(cellstr(num2str(cellNum'))'); % convert to cell array of strings cellNum = {'1' '2' '3'}
        tc_temp = flipud(tc_temp); % flip the data up-down for reverse plot order along y axis
end

x_end = size(tc_temp,2);

clear tch;
for iC = 1:size(tc_temp,1)      
    
    subaxis(size(tc_temp,1),1,iC,'SpacingVert',0);
    
    switch cfg.mode
        case 'area'
            tch(iC) = area(tc_temp(iC,:)); set(tch(iC),'FaceColor',cfg.color,'EdgeColor',[1 1 1]); 
            if cfg.alpha > 0 % attempt at transparency doesn't work
                patchobjs = findobj(get(tch(iC),'children'), '-property', 'AlphaData');
                set(patchobjs, 'AlphaData', cfg.alpha);
            end
        case 'line'
            tch(iC) = plot(tc_temp(iC,:),'Color',cfg.color,'LineWidth',1);
    end
    hold on; 
        
    %ylims = get(gca,'Ylim');
    %plot([TC.peak_idx(iC) TC.peak_idx(iC)],[ylims(1) ylims(2)],'r:');
    
    set(gca, 'XTick', [],'YTick',[],'XLim',[1 x_end],'LineWidth',0.5,'XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5],'XDir',cfg.XDir); 
    box off;
    
    yl = ylabel(cellNum(iC)); 
    set(yl,'Rotation',0,'Fontsize',8);
    %if iC == 1; title(ENC_data(iT).trial_type); end
    %axis off;
end

set(gca,'XTick',[1 cfg.cp x_end],'XTicklabel',[1 cfg.cp*cfg.binsize x_end*cfg.binsize],'Ticklength', [0 0],'Fontsize',8,'YTickLabel',cellNum);


