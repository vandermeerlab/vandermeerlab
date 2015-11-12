% POSTER_Plot_CoOccur

% this takes in the output from the SfN rushcode script
% POSTER_CollectCorrSequences

% generates a multiplot of all sequences data. Choose whether you want 
% matched or non matche dfield data plotted (data must exist already). 
%
% aacarey Oct 2015

clear

% WHAT DO YOU WANT THIS SCRIPT TO DO?

% some appearance things
cfg.colormode = 'inventory3';
FontSize = 8;

cfg.useMatched = 1; % use matched fields data

% where to save the figure?
cfg.output_fd = 'E:\Documents\TmazePaper\visuals';

cfg.showAllRatsText = 1; % do you want to show the text "all rats" on the combined data figures?

%%
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd('E:\Documents\TmazePaper\data')
if cfg. useMatched
    pre = load('CS_PRE_matched_sequences_rush.mat');
    task = load('CS_TASK_matched_sequences_rush.mat');
    post = load('CS_POST_matched_sequences_rush.mat');
    cfg.output_fn = 'seqMulti_matched_rush';
else
    pre = load('CS_PRE_sequences_rush.mat');
    task = load('CS_TASK_sequences_rush.mat');
    post = load('CS_POST_sequences_rush.mat');
    cfg.output_fn = 'seqMulti_rush';
end

ylab = {'Number of significant'; 'sequences'};
ylimsall = [0 110];
yticksall = ylimsall(1):10:ylimsall(2);
ylimssing = [0 80];
ytickssing = ylimssing(1):20:ylimssing(2);
cfg.output_fn = 'seqMulti_matched_rush';



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
cd(cfg.output_fd); 
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-dpdf','-r300',cfg.output_fn);
print(gcf,'-depsc','-r300',cfg.output_fn);
cd(originalFolder)