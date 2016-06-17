% ALL_Plot_DeqSeq
%
% plot output from ALL_Collect_DecSeq

% some appearance things
cfg.colormode = 'inventory3';
FontSize = 8;
cfg.input_fd = 'D:\projects\AlyssaTmaze\resultsFiles\';
cfg.output_fd = 'D:\projects\AlyssaTmaze\resultsFiles\viz';
cfg.showAllRatsText = 1; % do you want to show the text "all rats" on the combined data figures?
cfg.writeOutput = 1;
cfg.outbasefn = 'DecSeqMulti'; % base filename for figure output
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;


%% F1: raw sequence counts
%cfg.ylim = [600 300]; cfg.ylimtick = [150 75]; % ALL - overall and single lims
cfg.ylim = [300 150]; cfg.ylimtick = [75 37.5]; % ALL - overall and single lims

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_counts');

pre = load('DecSeq_prerecord_all_out');
task = load('DecSeq_taskrest_all_out');
post = load('DecSeq_postrecord_all_out');

ylab = {'Number of significant'; 'sequences'};
ylimsall = [0 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [0 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE
d = [pre.data.(rats{iRat}).food_left pre.data.(rats{iRat}).food_right pre.data.(rats{iRat}).water_left pre.data.(rats{iRat}).water_right];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',yticksall)
xlabel('  food                  water','FontSize',FontSize)
ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
title('PRERECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[3 4 9 10]) % for all rats TASKREST
d = [task.data.(rats{iRat}).food_left task.data.(rats{iRat}).food_right task.data.(rats{iRat}).water_left task.data.(rats{iRat}).water_right];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('TASK')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[5 6 11 12]) % for all rats POST
d = [post.data.(rats{iRat}).food_left post.data.(rats{iRat}).food_right post.data.(rats{iRat}).water_left post.data.(rats{iRat}).water_right];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('POSTRECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.65; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    d = [pre.data.(rats{iRat}).food_left pre.data.(rats{iRat}).food_right pre.data.(rats{iRat}).water_left pre.data.(rats{iRat}).water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    if position(1) == 13 || position(1) == 19
        yticks = ytickssing;
        ylabel(ylab,'FontSize',FontSize)
    else
        yticks = [];
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',yticks)
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(2)) % for indiv rat TASK
    d = [task.data.(rats{iRat}).food_left task.data.(rats{iRat}).food_right task.data.(rats{iRat}).water_left task.data.(rats{iRat}).water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(3)) % for indiv rat POST
    d = [post.data.(rats{iRat}).food_left post.data.(rats{iRat}).food_right post.data.(rats{iRat}).water_left post.data.(rats{iRat}).water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-dpdf','-r300',cfg.output_fn);
    print(gcf,'-depsc','-r300',cfg.output_fn);
    cd(originalFolder)
end

%% F2: proportional sequence counts
cfg.ylim = [1 1]; cfg.ylimtick = [0.25 0.25]; % overall and single lims

colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_props');

ylab = {'Proportion of significant'; 'sequences'};
ylimsall = [0 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [0 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE
d = [pre.data.(rats{iRat}).food_leftN pre.data.(rats{iRat}).food_rightN pre.data.(rats{iRat}).water_leftN pre.data.(rats{iRat}).water_rightN];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',yticksall)
xlabel('  food                  water','FontSize',FontSize)
ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
title('PRERECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[3 4 9 10]) % for all rats TASKREST
d = [task.data.(rats{iRat}).food_leftN task.data.(rats{iRat}).food_rightN task.data.(rats{iRat}).water_leftN task.data.(rats{iRat}).water_rightN];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('TASK')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[5 6 11 12]) % for all rats POST
d = [post.data.(rats{iRat}).food_leftN post.data.(rats{iRat}).food_rightN post.data.(rats{iRat}).water_leftN post.data.(rats{iRat}).water_rightN];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('POSTRECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.65; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    d = [pre.data.(rats{iRat}).food_leftN pre.data.(rats{iRat}).food_rightN pre.data.(rats{iRat}).water_leftN pre.data.(rats{iRat}).water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    if position(1) == 13 || position(1) == 19
        yticks = ytickssing;
        ylabel(ylab,'FontSize',FontSize)
    else
        yticks = [];
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',yticks)
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(2)) % for indiv rat TASK
    d = [task.data.(rats{iRat}).food_leftN task.data.(rats{iRat}).food_rightN task.data.(rats{iRat}).water_leftN task.data.(rats{iRat}).water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(3)) % for indiv rat POST
    d = [post.data.(rats{iRat}).food_leftN post.data.(rats{iRat}).food_rightN post.data.(rats{iRat}).water_leftN post.data.(rats{iRat}).water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto');  drawnow;
    cd(cfg.output_fd);
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-dpdf','-r300',cfg.output_fn);
    print(gcf,'-depsc','-r300',cfg.output_fn);
    cd(originalFolder)
end

%% F3: relative differences
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_props_rel');

cfg.xlim = [0.5 0.8]; cfg.xlimtick = [0.5 0.8]; % overall and single lims

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).w}; % {foodColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE

food_diff = pre.data.(rats{iRat}).food_rightN-pre.data.(rats{iRat}).food_leftN;
water_diff = pre.data.(rats{iRat}).water_rightN-pre.data.(rats{iRat}).water_leftN;

h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
hold on
h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(1) cfg.xlim(1)],'XTick',-cfg.xlim(1):cfg.xlimtick(1)/2:cfg.xlim(1), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(1)),-cfg.xlim(1)/2,0,cfg.xlim(1)/2,sprintf('W %.1f',cfg.xlim(1))}, ...
    'YTick',1:2,'YTickLabel',{'WaterR','FoodR'},'FontSize',FontSize,'YLim',[0.5 2.5]);

title('PRERECORD')
box off
set(gca,'Layer','top')

if cfg.showAllRatsText
txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[3 4 9 10]) % for all rats TASKREST
food_diff = task.data.(rats{iRat}).food_rightN-task.data.(rats{iRat}).food_leftN;
water_diff = task.data.(rats{iRat}).water_rightN-task.data.(rats{iRat}).water_leftN;

h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
hold on
h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(1) cfg.xlim(1)],'XTick',-cfg.xlim(1):cfg.xlimtick(1)/2:cfg.xlim(1), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(1)),-cfg.xlim(1)/2,0,cfg.xlim(1)/2,sprintf('W %.1f',cfg.xlim(1))}, ...
    'YTick',1:2,'YTickLabel',{'WaterR','RoodR'},'FontSize',FontSize,'YLim',[0.5 2.5]);

title('TASK')
box off
set(gca,'Layer','top')

if cfg.showAllRatsText
txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[5 6 11 12]) % for all rats POST

food_diff = post.data.(rats{iRat}).food_rightN-post.data.(rats{iRat}).food_leftN;
water_diff = post.data.(rats{iRat}).water_rightN-post.data.(rats{iRat}).water_leftN;

h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
hold on
h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(1) cfg.xlim(1)],'XTick',-cfg.xlim(1):cfg.xlimtick(1)/2:cfg.xlim(1), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(1)),-cfg.xlim(1)/2,0,cfg.xlim(1)/2,sprintf('W %.1f',cfg.xlim(1))}, ...
    'YTick',1:2,'YTickLabel',{'WaterR','FoodR'},'FontSize',FontSize,'YLim',[0.5 2.5]);

title('POST')
box off
set(gca,'Layer','top')

if cfg.showAllRatsText
txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.2; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).w}; % {foodColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    food_diff = pre.data.(rats{iRat}).food_rightN-pre.data.(rats{iRat}).food_leftN;
    water_diff = pre.data.(rats{iRat}).water_rightN-pre.data.(rats{iRat}).water_leftN;
    
    h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
    hold on
    h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(2) cfg.xlim(2)],'XTick',-cfg.xlim(2):cfg.xlimtick(2)/2:cfg.xlim(2), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(2)),'',0,'',sprintf('W %.1f',cfg.xlim(2))}, ...
    'YTick',1:2,'YTickLabel',{'Wr','Fr'},'FontSize',FontSize,'YLim',[0.5 2.5]);
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(2)) % for indiv rat TASK
    
    food_diff = task.data.(rats{iRat}).food_rightN-task.data.(rats{iRat}).food_leftN;
    water_diff = task.data.(rats{iRat}).water_rightN-task.data.(rats{iRat}).water_leftN;
    
    h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
    hold on
    h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(2) cfg.xlim(2)],'XTick',-cfg.xlim(2):cfg.xlimtick(2)/2:cfg.xlim(2), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(2)),'',0,'',sprintf('W %.1f',cfg.xlim(2))}, ...
    'YTick',1:2,'YTickLabel',{'Wr','Fr'},'FontSize',FontSize,'YLim',[0.5 2.5]);
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(3)) % for indiv rat POST
    
    food_diff = post.data.(rats{iRat}).food_rightN-post.data.(rats{iRat}).food_leftN;
    water_diff = post.data.(rats{iRat}).water_rightN-post.data.(rats{iRat}).water_leftN;
    
    h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
    hold on
    h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(2) cfg.xlim(2)],'XTick',-cfg.xlim(2):cfg.xlimtick(2)/2:cfg.xlim(2), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(2)),'',0,'',sprintf('W %.1f',cfg.xlim(2))}, ...
    'YTick',1:2,'YTickLabel',{'Wr','Fr'},'FontSize',FontSize,'YLim',[0.5 2.5]);
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    cd(cfg.output_fd);
    maximize; set(gcf,'PaperPositionMode','auto');  drawnow;
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-dpdf','-r300',cfg.output_fn);
    print(gcf,'-depsc','-r300',cfg.output_fn);
    cd(originalFolder)
end


%% F4: decoding error
cfg.ylim = [20 30]; cfg.ylimtick = [5 10]; % overall and single lims

colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_decAcc');

pre = load('DecSeq_prerecord_all_out');

ylab = {'decoding error'};

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE

this_data = pre.data.(rats{iRat}).ALL_sig_seq;

food_left_idx = find(this_data.arm == 1 & this_data.type == 1);
d(1) = nanmean(this_data.decErr(food_left_idx));
food_right_idx = find(this_data.arm == 2 & this_data.type == 1);
d(2) = nanmean(this_data.decErr(food_right_idx));
water_left_idx = find(this_data.arm == 1 & this_data.type == 2);
d(3) = nanmean(this_data.decErr(water_left_idx));
water_right_idx = find(this_data.arm == 2 & this_data.type == 2);
d(4) = nanmean(this_data.decErr(water_right_idx));

for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',[0 cfg.ylim(1)],'YTick',0:cfg.ylimtick(1):cfg.ylim(1));
xlabel('  food                  water','FontSize',FontSize)
ylabel(ylab,'FontSize',FontSize)
title('DecAcc')
box off

set(gca,'Layer','top')
if cfg.showAllRatsText
    txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.2; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
      
    subplot(4,6,position(1)) % for indiv rat PRE
    
    
    this_data = pre.data.(rats{iRat}).ALL_sig_seq;
    
    food_left_idx = find(this_data.arm == 1 & this_data.type == 1);
    d(1) = nanmean(this_data.decErr(food_left_idx));
    food_right_idx = find(this_data.arm == 2 & this_data.type == 1);
    d(2) = nanmean(this_data.decErr(food_right_idx));
    water_left_idx = find(this_data.arm == 1 & this_data.type == 2);
    d(3) = nanmean(this_data.decErr(water_left_idx));
    water_right_idx = find(this_data.arm == 2 & this_data.type == 2);
    d(4) = nanmean(this_data.decErr(water_right_idx));
    
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',[0 cfg.ylim(2)],'YTick',0:cfg.ylimtick(2):cfg.ylim(2));
    xlabel('  food                  water','FontSize',FontSize)
    ylabel(ylab,'FontSize',FontSize)
    box off
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    cd(cfg.output_fd);
    maximize; set(gcf,'PaperPositionMode','auto');  drawnow;
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-dpdf','-r300',cfg.output_fn);
    print(gcf,'-depsc','-r300',cfg.output_fn);
    cd(originalFolder)
end

%% F5: arranged by first choice
cfg.ylim = [1 1]; cfg.ylimtick = [0.25 0.25]; % overall and single lims

colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_firstchoice');

ylab = {'Proportion of significant'; 'sequences'};
ylimsall = [0 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [0 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE

this_data = pre.data.(rats{iRat}).ALL_sig_seq;

choiceleft_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 1);
choiceleft_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 1);

choiceright_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 2);
choiceright_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 2);


xbin = [1 2 4 5];
d = [nanmean(choiceleft_left) nanmean(choiceleft_right) nanmean(choiceright_left) nanmean(choiceright_right)];
cl = 'rrbb';

for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',yticksall)
xlabel('  left                  right','FontSize',FontSize)
ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
title('PRERECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.2; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    this_data = pre.data.(rats{iRat}).ALL_sig_seq;
    
    choiceleft_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 1);
    choiceleft_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 1);
    
    choiceright_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 2);
    choiceright_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 2);
    
    
    xbin = [1 2 4 5];
    d = [nanmean(choiceleft_left) nanmean(choiceleft_right) nanmean(choiceright_left) nanmean(choiceright_right)];
    cl = 'rrbb';
    
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',yticksall)
    xlabel('  left                  right','FontSize',FontSize)
    ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    cd(cfg.output_fd);
    maximize; set(gcf,'PaperPositionMode','auto');  drawnow;
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-dpdf','-r300',cfg.output_fn);
    print(gcf,'-depsc','-r300',cfg.output_fn);
    cd(originalFolder)
end