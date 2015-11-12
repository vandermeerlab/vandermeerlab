%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                        The Script that Rhymes                       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Tmaze replays beside tuning curves.

% aacarey Sept 2015


clear

% what to call the output figure and where to save it
cfg.output_fn = 'F3test'; % it saves as a .png
cfg.output_fd = 'E:\Documents\TmazePaper\visuals';

cfg.FontSize = 6;
cfg.HorzLoc = 0; % in normalized units, the percentage of total x axis breadth that the rat's ID will sit at

cfg.colormode = 'inventory3';

col = TmazeColors(cfg.colormode);

% PLOT R042
%% load some stuff
cd('E:\data\promoted\R042\R042-2013-08-18')
t = 5042.7;

iWasHere = pwd; % remember where I started

LoadExpKeys % loads variable 'ExpKeys'
LoadCandidates % loads variable 'evt'

% grab the pc from \files
cd([pwd,'\files'])
load(FindFile('*-armPC.mat')) % loads variable named 'armPC'
cd(iWasHere)

cfg_temp = [];
cfg_temp.load_questionable_cells = 1;
S = LoadSpikes(cfg_temp);

cfg_temp = [];
cfg_temp.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(cfg_temp);

%%

gap = [0.01,0.01]; % for subtightplot [vertical horizontal] spacing b/t subs
marg_h = [];
marg_w = [0.09 0.06];
% set up multiraster configs
mr1 = [];
mr1.lfp = csc;
mr1.lfpHeight = 10;
mr1.evt = evt;
mr1.ivColor = 'k';
mr1.openNewFig = 0;
mr1.SpikeHeight = 0.5;
mr1.LineWidth = 0.5;
mr1.spkColor = col.all.f;
mr1.setAxes = 'on';

mr2 = mr1;
mr2.spkColor = col.all.w;

% set up tuning curve configs

cfg_tc1 = [];
%cfg_tc.arrange = 'stack'; % not used for heat mode
cfg_tc1.mode = 'heat';
cfg_tc1.select = 1;
cfg_tc1.LineWidth = 1;
cfg_tc1.YOrd = 'normal';
cfg_tc1.color = 'hot';
cfg_tc1.EdgeColor = 'k';
cfg_tc1.ymin = -9;
cfg_tc1.openNewFig = 0;

cfg_tc2 = cfg_tc1;
cfg_tc2.ymin = -10;
% plot some stuff

figure('KeyPressFcn',@navigate)


ax1 = subtightplot(15,4,[1 5],gap,marg_h,marg_w);

MultiRaster(mr1,armPC.L.S);
%set(ax1,'Color','none')
set(ax1,'YTickLabel',[],'XTickLabel',[])
set(ax1,'ticklength',[0 0],'ycolor','w','xcolor','w')
ylabel('Left','Color','k','FontSize',cfg.FontSize)
whylimbs = get(gca,'YLim');
set(gca,'YLim',[whylimbs(1) whylimbs(2)+mr1.SpikeHeight])

txt = 'R042: right sequence, water restricted';
    text(cfg.HorzLoc,1,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','left',...
        'FontSize',cfg.FontSize)

ax4 = subtightplot(15,4,[10 14],gap,marg_h,marg_w);
cfg_tc1.openNewFig = 0;
cfg_tc1.ymin = whylimbs(1);
PlotTC2(cfg_tc2,armPC.R.tc);
set(ax4,'ycolor','k','xcolor','w','yaxislocation','right','yticklabel',0:length(armPC.R.S.t):length(armPC.R.S.t),'YTick',0:length(armPC.R.S.t):length(armPC.R.S.t),'XTickLabel',[],'FontSize',cfg.FontSize)
%set(ax4,'ycolor','w','xcolor','w')
text(1.07,0.77,'cell #','Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.FontSize,'rot',90)
    mafak = get(gca,'YLim');
choicepoint = size(armPC.R.tc.tc,2)-42;
plot([choicepoint choicepoint],mafak,'w','LineWidth',2,'LineStyle',':')
ax4b = axes('Position', get(ax4, 'Position'),'Color','none');
set(ax4b,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')

%xlims = get(ax4,'XLim');
%xPath = size(armPC.R.tc.tc,2)*3;
%xConv = xPath/(xlims(2)-xlims(1));


ax2 = subtightplot(15,4,[2 6],gap,marg_h,marg_w);

PlotTC2(cfg_tc1,armPC.L.tc);
set(ax2,'ycolor','k','xcolor','w','yaxislocation','right','yticklabel',0:length(armPC.L.S.t):length(armPC.L.S.t),'YTick',0:length(armPC.L.S.t):length(armPC.L.S.t),'XTickLabel',[],'FontSize',cfg.FontSize)
text(1.07,0.77,'cell #','Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.FontSize,'rot',90)
mafak = get(gca,'YLim');
choicepoint = size(armPC.R.tc.tc,2)-42;
plot([choicepoint choicepoint],mafak,'w','LineWidth',2,'LineStyle',':')
ax2b = axes('Position', get(ax2, 'Position'),'Color','none');
set(ax2b,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')

ax3 = subtightplot(15,4,[9 13],gap,marg_h,marg_w);
%hold on;
MultiRaster(mr2,armPC.R.S);
%set(ax3,'Color','none')
set(ax3,'YTickLabel',[],'XTickLabel',[])
set(ax3,'ticklength',[0 0],'ycolor','w','xcolor','w')
ylabel('Right','Color','k','FontSize',cfg.FontSize)
whylimbs = get(gca,'YLim');
set(gca,'YLim',[whylimbs(1) whylimbs(2)+mr1.SpikeHeight])


linkaxes([ax1 ax3],'x'); % make it so i can navigate through both raster plots at the same time

%%

set(ax1,'XLim',[t-0.5 t+0.5]);

%% ADD SCALE BARS

ax5 = subtightplot(15,4,17,gap,marg_h,marg_w);
% make a time scalebar like a non-pro
yhorz = [1 1];
x = [9 10];
yvert = [0.96 1.04];
plot(x,yhorz,'Color','k','LineWidth',1); hold on;
plot([x(1) x(1)],yvert,'Color','k','LineWidth',1)
plot([x(2) x(2)],yvert,'Color','k','LineWidth',1)
set(gca,'XLim',[0 10],'YLim',[0.5 yvert(2)])% ,'Position',[0 0.9 0.9 0.1])
% annotate that thing
xSpot = 8.7; ySpot = 0.75;
text(xSpot,ySpot,'100 ms','FontSize',cfg.FontSize)
box off
axis off

ax6 = subtightplot(15,4,18,gap,marg_h,marg_w);
xlims = get(ax4,'XLim');
%posBinSize = 3; % 3, how many cm per bin (when position data was standardized)
yhorz = [1 1];
x = [xlims(2)-10 xlims(2)];
plot(x,yhorz,'Color','k','LineWidth',1); hold on;
plot([x(1) x(1)],yvert,'Color','k','LineWidth',1)
plot([x(2) x(2)],yvert,'Color','k','LineWidth',1)
set(gca,'XLim',xlims,'YLim',[0.5 yvert(2)])% ,'Position',[0 0.9 0.9 0.1])
% annotate that thing
xSpot = xlims(2)*0.87; ySpot = 0.75;
text(xSpot,ySpot,'30 cm','FontSize',cfg.FontSize)
box off
axis off

%% PLOT R050
%% load some stuff
cd('E:\data\promoted\R050\R050-2014-04-02')
t = 4563.2;

iWasHere = pwd; % remember where I started

LoadExpKeys % loads variable 'ExpKeys'
LoadCandidates % loads variable 'evt'

% grab the pc from \files
cd([pwd,'\files'])
load(FindFile('*-armPC.mat')) % loads variable named 'armPC'
cd(iWasHere)

cfg_temp = [];
cfg_temp.load_questionable_cells = 1;
S = LoadSpikes(cfg_temp);

cfg_temp = [];
cfg_temp.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(cfg_temp);

%%

gap = [0.01,0.01]; % for subtightplot [vertical horizontal] spacing b/t subs
marg_h = [];
marg_w = [0.09 0.06];
% set up multiraster configs
mr1 = [];
mr1.lfp = csc;
mr1.lfpHeight = 10;
mr1.lfpMax = 6;
mr1.evt = evt;
mr1.ivColor = 'k';
mr1.openNewFig = 0;
mr1.SpikeHeight = 0.5;
mr1.LineWidth = 0.5;
mr1.spkColor = col.all.f;
mr1.setAxes = 'on';

mr2 = mr1;
mr2.spkColor = col.all.w;

% set up tuning curve configs

cfg_tc1 = [];
%cfg_tc.arrange = 'stack'; % not used for heat mode
cfg_tc1.mode = 'heat';
cfg_tc1.select = 1;
cfg_tc1.LineWidth = 1;
cfg_tc1.YOrd = 'normal';
cfg_tc1.color = 'hot';
cfg_tc1.EdgeColor = 'k';
cfg_tc1.ymin = -9;
cfg_tc1.openNewFig = 0;

cfg_tc2 = cfg_tc1;
cfg_tc2.ymin = -10;
% plot some stuff


ax1 = subtightplot(15,4,[21 25],gap,marg_h,marg_w);

MultiRaster(mr1,armPC.L.S);
%set(ax1,'Color','none')
set(ax1,'YTickLabel',[],'XTickLabel',[])
set(ax1,'ticklength',[0 0],'ycolor','w','xcolor','w')
ylabel('Left','Color','k','FontSize',cfg.FontSize)
whylimbs = get(gca,'YLim');
set(gca,'YLim',[whylimbs(1) whylimbs(2)+mr1.SpikeHeight])

txt = 'R050: right sequence, food restricted';
    text(cfg.HorzLoc,1,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','left',...
        'FontSize',cfg.FontSize)

ax4 = subtightplot(15,4,[30 34],gap,marg_h,marg_w);
cfg_tc1.openNewFig = 0;
cfg_tc1.ymin = whylimbs(1);
PlotTC2(cfg_tc2,armPC.R.tc);
mafak = get(gca,'YLim');
choicepoint = size(armPC.R.tc.tc,2)-61;
plot([choicepoint choicepoint],mafak,'w','LineWidth',2,'LineStyle',':')
set(ax4,'ycolor','k','xcolor','w','yaxislocation','right','yticklabel',0:length(armPC.R.S.t):length(armPC.R.S.t),'YTick',0:length(armPC.R.S.t):length(armPC.R.S.t),'XTickLabel',[],'FontSize',cfg.FontSize)
%set(ax4,'ycolor','w','xcolor','w')
text(1.07,0.77,'cell #','Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.FontSize,'rot',90)

ax4b = axes('Position', get(ax4, 'Position'),'Color','none');
set(ax4b,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')

%xlims = get(ax4,'XLim');
%xPath = size(armPC.R.tc.tc,2)*3;
%xConv = xPath/(xlims(2)-xlims(1));


ax2 = subtightplot(15,4,[22 26],gap,marg_h,marg_w);

PlotTC2(cfg_tc1,armPC.L.tc);
set(ax2,'ycolor','k','xcolor','w','yaxislocation','right','yticklabel',0:length(armPC.L.S.t):length(armPC.L.S.t),'YTick',0:length(armPC.L.S.t):length(armPC.L.S.t),'XTickLabel',[],'FontSize',cfg.FontSize)
text(1.07,0.77,'cell #','Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.FontSize,'rot',90)
mafak = get(gca,'YLim');
choicepoint = size(armPC.R.tc.tc,2)-61;
plot([choicepoint choicepoint],mafak,'w','LineWidth',2,'LineStyle',':')
ax2b = axes('Position', get(ax2, 'Position'),'Color','none');
set(ax2b,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')


ax3 = subtightplot(15,4,[29 33],gap,marg_h,marg_w);
%hold on;
MultiRaster(mr2,armPC.R.S);
%set(ax3,'Color','none')
set(ax3,'YTickLabel',[],'XTickLabel',[])
set(ax3,'ticklength',[0 0],'ycolor','w','xcolor','w')
ylabel('Right','Color','k','FontSize',cfg.FontSize)
whylimbs = get(gca,'YLim');
set(gca,'YLim',[whylimbs(1) whylimbs(2)+mr1.SpikeHeight])


linkaxes([ax1 ax3],'x'); % make it so i can navigate through both raster plots at the same time

%%

set(ax1,'XLim',[t-0.5 t+0.5]);
%% ADD SCALE BARS

ax5 = subtightplot(15,4,37,gap,marg_h,marg_w);
% make a time scalebar like a non-pro
yhorz = [1 1];
x = [9 10];
yvert = [0.96 1.04];
plot(x,yhorz,'Color','k','LineWidth',1); hold on;
plot([x(1) x(1)],yvert,'Color','k','LineWidth',1)
plot([x(2) x(2)],yvert,'Color','k','LineWidth',1)
set(gca,'XLim',[0 10],'YLim',[0.5 yvert(2)])% ,'Position',[0 0.9 0.9 0.1])
% annotate that thing
xSpot = 8.7; ySpot = 0.75;
text(xSpot,ySpot,'100 ms','FontSize',cfg.FontSize)
box off
axis off

ax6 = subtightplot(15,4,38,gap,marg_h,marg_w);
% make a position scalebar 
xlims = get(ax4,'XLim');
%posBinSize = 3; % 3, how many cm per bin (when position data was standardized)
yhorz = [1 1];
x = [xlims(2)-10 xlims(2)];
plot(x,yhorz,'Color','k','LineWidth',1); hold on;
plot([x(1) x(1)],yvert,'Color','k','LineWidth',1)
plot([x(2) x(2)],yvert,'Color','k','LineWidth',1)
set(gca,'XLim',xlims,'YLim',[0.5 yvert(2)])% ,'Position',[0 0.9 0.9 0.1])
% annotate that thing
xSpot = xlims(2)*0.89; ySpot = 0.75;
text(xSpot,ySpot,'30 cm','FontSize',cfg.FontSize)
box off
axis off
%% PLOT R064
%% load some stuff
cd('E:\data\promoted\R064\R064-2015-04-23')
t = 3191.4;

iWasHere = pwd; % remember where I started

LoadExpKeys % loads variable 'ExpKeys'
LoadCandidates % loads variable 'evt'

% grab the pc from \files
cd([pwd,'\files'])
load(FindFile('*-armPC.mat')) % loads variable named 'armPC'
cd(iWasHere)

cfg_temp = [];
cfg_temp.load_questionable_cells = 1;
S = LoadSpikes(cfg_temp);

cfg_temp = [];
cfg_temp.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(cfg_temp);

%%

gap = [0.01,0.01]; % for subtightplot [vertical horizontal] spacing b/t subs
marg_h = [];
marg_w = [0.09 0.06];
% set up multiraster configs
mr1 = [];
mr1.lfp = csc;
mr1.lfpHeight = 10;
mr1.lfpMax = 5;
mr1.evt = evt;
mr1.ivColor = 'k';
mr1.openNewFig = 0;
mr1.SpikeHeight = 0.5;
mr1.LineWidth = 0.5;
mr1.spkColor = col.all.f;
mr1.setAxes = 'on';

mr2 = mr1;
mr2.spkColor = col.all.w;

% set up tuning curve configs

cfg_tc1 = [];
%cfg_tc.arrange = 'stack'; % not used for heat mode
cfg_tc1.mode = 'heat';
cfg_tc1.select = 1;
cfg_tc1.LineWidth = 1;
cfg_tc1.YOrd = 'normal';
cfg_tc1.color = 'hot';
cfg_tc1.EdgeColor = 'k';
cfg_tc1.ymin = -9;
cfg_tc1.openNewFig = 0;

cfg_tc2 = cfg_tc1;
cfg_tc2.ymin = -10;
% plot some stuff


ax1 = subtightplot(15,4,[41 45],gap,marg_h,marg_w);

MultiRaster(mr1,armPC.L.S);
%set(ax1,'Color','none')
set(ax1,'YTickLabel',[],'XTickLabel',[])
set(ax1,'ticklength',[0 0],'ycolor','w','xcolor','w')
ylabel('Left','Color','k','FontSize',cfg.FontSize)
whylimbs = get(gca,'YLim');
set(gca,'YLim',[whylimbs(1) whylimbs(2)+mr1.SpikeHeight])

txt = 'R064: left sequence, water restricted';
    text(cfg.HorzLoc,1,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','left',...
        'FontSize',cfg.FontSize)


ax4 = subtightplot(15,4,[50 54],gap,marg_h,marg_w);
cfg_tc1.openNewFig = 0;
cfg_tc1.ymin = whylimbs(1);
PlotTC2(cfg_tc2,armPC.R.tc);
set(ax4,'ycolor','k','xcolor','w','yaxislocation','right','yticklabel',0:length(armPC.R.S.t):length(armPC.R.S.t),'YTick',0:length(armPC.R.S.t):length(armPC.R.S.t),'XTickLabel',[],'FontSize',cfg.FontSize)
%set(ax4,'ycolor','w','xcolor','w')
text(1.07,0.77,'cell #','Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.FontSize,'rot',90)
    mafak = get(gca,'YLim');
choicepoint = size(armPC.R.tc.tc,2)-61;
plot([choicepoint choicepoint],mafak,'w','LineWidth',2,'LineStyle',':')
ax4b = axes('Position', get(ax4, 'Position'),'Color','none');
set(ax4b,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')

%xlims = get(ax4,'XLim');
%xPath = size(armPC.R.tc.tc,2)*3;
%xConv = xPath/(xlims(2)-xlims(1));


ax2 = subtightplot(15,4,[42 46],gap,marg_h,marg_w);

PlotTC2(cfg_tc1,armPC.L.tc);

set(ax2,'ycolor','k','xcolor','w','yaxislocation','right','yticklabel',0:length(armPC.L.S.t):length(armPC.L.S.t),'YTick',0:length(armPC.L.S.t):length(armPC.L.S.t),'XTickLabel',[],'FontSize',cfg.FontSize)
text(1.07,0.77,'cell #','Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',cfg.FontSize,'rot',90)
mafak = get(gca,'YLim');
choicepoint = size(armPC.R.tc.tc,2)-61;
plot([choicepoint choicepoint],mafak,'w','LineWidth',2,'LineStyle',':')
ax2b = axes('Position', get(ax2, 'Position'),'Color','none');
set(ax2b,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')


ax3 = subtightplot(15,4,[49 53],gap,marg_h,marg_w);
%hold on;
MultiRaster(mr2,armPC.R.S);
%set(ax3,'Color','none')
set(ax3,'YTickLabel',[],'XTickLabel',[])
set(ax3,'ticklength',[0 0],'ycolor','w','xcolor','w')
ylabel('Right','Color','k','FontSize',cfg.FontSize)
whylimbs = get(gca,'YLim');
set(gca,'YLim',[whylimbs(1) whylimbs(2)+mr1.SpikeHeight])


linkaxes([ax1 ax3],'x'); % make it so i can navigate through both raster plots at the same time

%%

set(ax1,'XLim',[t-0.5 t+0.5]);

%% ADD SCALE BARS

ax5 = subtightplot(15,4,57,gap,marg_h,marg_w);
% make a time scalebar like a non-pro
yhorz = [1 1];
x = [9 10];
yvert = [0.96 1.04];
plot(x,yhorz,'Color','k','LineWidth',1); hold on;
plot([x(1) x(1)],yvert,'Color','k','LineWidth',1)
plot([x(2) x(2)],yvert,'Color','k','LineWidth',1)
set(gca,'XLim',[0 10],'YLim',[0.5 yvert(2)])% ,'Position',[0 0.9 0.9 0.1])
% annotate that thing
xSpot = 8.7; ySpot = 0.75;
text(xSpot,ySpot,'100 ms','FontSize',cfg.FontSize)
box off
axis off

ax6 = subtightplot(15,4,58,gap,marg_h,marg_w);
% make a position scalebar 
xlims = get(ax4,'XLim');
%posBinSize = 3; % 3, how many cm per bin (when position data was standardized)
yhorz = [1 1];
x = [xlims(2)-10 xlims(2)];
plot(x,yhorz,'Color','k','LineWidth',1); hold on;
plot([x(1) x(1)],yvert,'Color','k','LineWidth',1)
plot([x(2) x(2)],yvert,'Color','k','LineWidth',1)
set(gca,'XLim',xlims,'YLim',[0.5 yvert(2)])% ,'Position',[0 0.9 0.9 0.1])
% annotate that thing
xSpot = xlims(2)*0.89; ySpot = 0.75;
text(xSpot,ySpot,'30 cm','FontSize',cfg.FontSize)
box off
axis off

%%
set(gcf,'Color','w')
set(gcf,'InvertHardCopy','off')
cd(cfg.output_fd)
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-depsc','-r300',cfg.output_fn);
print(gcf,'-dpdf','-r300',cfg.output_fn);

cd(iWasHere)
