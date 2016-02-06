% POSTER_Plot_CoOccur

% this takes in the output from the SfN rushcode script POSTER_Collect_CoOccur

% generates a multiplot of all co-occurrence data for both p0 and p4 
%

clear

% WHAT DO YOU WANT THIS SCRIPT TO DO?

% what to plot?
cfg.whichP = 'p4'; % 'p0' or 'p4' (single cell or zscore coactivity)

% some appearance things
cfg.colormode = 'inventory3';
FontSize = 8;

% where to save the figure?
cfg.input_fd = 'E:\Documents\TmazePaper\data';
cfg.output_fd = 'E:\Documents\TmazePaper\visuals';

%%
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)

pre = load('CC_PRE_cooccurrence_rush.mat');
task = load('CC_TASK_cooccurrence_rush.mat');
post = load('CC_POST_cooccurrence_rush.mat');

switch cfg.whichP
    case 'p4'
        ylab = {'SWR coactivation'; 'Z-score'};
        ylimsall = [0 1.25];
        yticksall = ylimsall(1):0.25:ylimsall(2);
        ylimssing = [-1 2.5];
        ytickssing = ylimssing(1):ylimssing(2); 
        cfg.output_fn = 'p4Multi_rush';
        
    case 'p0'
        ylab = {'Proportion of';'SWRs active'};
        ylimsall = [0 0.35];
        yticksall = ylimsall(1):0.05:ylimsall(2);
        ylimssing = [0 0.35];
        ytickssing = ylimssing(1):0.1:ylimssing(2);
        cfg.output_fn = 'p0Multi_rush';
        
    otherwise
        error('Unrecognized cfg.whichP')
end

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];


% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    

figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE
d = [nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_food_left'])) nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_food_right'])) nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_water_left'])) nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_water_right']))];
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
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)

subplot(4,6,[3 4 9 10]) % for all rats TASKREST
d = [nanmean(task.data.(rats{iRat}).([cfg.whichP,'_food_left'])) nanmean(task.data.(rats{iRat}).([cfg.whichP,'_food_right'])) nanmean(task.data.(rats{iRat}).([cfg.whichP,'_water_left'])) nanmean(task.data.(rats{iRat}).([cfg.whichP,'_water_right']))];
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
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)

subplot(4,6,[5 6 11 12]) % for all rats POST
d = [nanmean(post.data.(rats{iRat}).([cfg.whichP,'_food_left'])) nanmean(post.data.(rats{iRat}).([cfg.whichP,'_food_right'])) nanmean(post.data.(rats{iRat}).([cfg.whichP,'_water_left'])) nanmean(post.data.(rats{iRat}).([cfg.whichP,'_water_right']))];
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
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.65; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    d = [nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_food_left'])) nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_food_right'])) nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_water_left'])) nanmean(pre.data.(rats{iRat}).([cfg.whichP,'_water_right']))];
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
    d = [nanmean(task.data.(rats{iRat}).([cfg.whichP,'_food_left'])) nanmean(task.data.(rats{iRat}).([cfg.whichP,'_food_right'])) nanmean(task.data.(rats{iRat}).([cfg.whichP,'_water_left'])) nanmean(task.data.(rats{iRat}).([cfg.whichP,'_water_right']))];
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
    d = [nanmean(post.data.(rats{iRat}).([cfg.whichP,'_food_left'])) nanmean(post.data.(rats{iRat}).([cfg.whichP,'_food_right'])) nanmean(post.data.(rats{iRat}).([cfg.whichP,'_water_left'])) nanmean(post.data.(rats{iRat}).([cfg.whichP,'_water_right']))];
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
cd(cfg.output_fd); 
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-depsc','-r300',cfg.output_fn);
print(gcf,'-dpdf','-r300',cfg.output_fn);
cd(originalFolder)