%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                      Plot Hypothesis Bar Graphs                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aacarey Oct 2015 for SfN poster

% choose some things about appearance
colors = TmazeColors('inventory3');
cfg.output_fn = 'hypothesis';
cfg.output_fd = 'E:\Documents\TmazePaper\visuals';

cfg.FontSize = 8;
cfg.LineWidth = 1;

xlims = [0 6];
ylims = [0 1];
%% pval branch thing nodes (will not make sense tomorrow)

x1 = 3; y1 = 0.95;
x2 = [1.5 3 4.5]; y2 = 0.9;
x3 = [1 2 4 5]; y3 = 0.85;
y4 = 0.8;

%% Do the thing

iWasHere = pwd;

xlab = '  food         water';

% numbers for "no effect of motivational state"
d1 = [0.55 0.45 0.52 0.48];
d2 = [0.4 0.6 0.4 0.6];

% numbers for "effect of motivational state"
d3 = [0.5 0.5 0.3 0.7];
d4 = [0.8 0.2 0.2 0.8];

col = {colors.all.f colors.all.f colors.all.w colors.all.w};
location = [1 2 4 5];

% plot some stuff

figure; hold on

% plot no effect ex 1
subplot(1,5,1)
for iBar = 1:length(d1)
    h1(iBar) = bar(location(iBar),d1(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XLim',xlims,'YLim',ylims,'YTick',[])
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',cfg.FontSize,'LineWidth',1)
xlabel(xlab)
ylabel('Neural representation')
pbaspect(gca,[1 1 1])
box off; 

text(1.2,1.1,'No effect of motivational state','Units','normalized','FontSize',cfg.FontSize,'HorizontalAlignment','center')

% plot no effect ex 2
subplot(1,5,2)
for iBar = 1:length(d2)
    h2(iBar) = bar(location(iBar),d2(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XLim',xlims,'YLim',ylims,'YTick',[])
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',cfg.FontSize,'LineWidth',1)
xlabel(xlab)
pbaspect(gca,[1 1 1])
box off; 


% plot effect ex 1
subplot(1,5,4)
for iBar = 1:length(d3)
    h3(iBar) = bar(location(iBar),d3(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XLim',xlims,'YLim',ylims,'YTick',[])
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',cfg.FontSize,'LineWidth',1)
xlabel(xlab)
ylabel('Neural representation')
pbaspect(gca,[1 1 1])
box off;

text(1.2,1.1,'Effect of motivational state','Units','normalized','FontSize',cfg.FontSize,'HorizontalAlignment','center')

% plot the pval branches
plot([x1 x1],[y1 y2],'k','LineWidth',cfg.LineWidth)
plot([x2(1) x2(3)],[y2 y2],'k','LineWidth',cfg.LineWidth)
plot([x2(1) x2(1)],[y2 y4],'k','LineWidth',cfg.LineWidth)
plot([x2(3) x2(3)],[y2 y3],'k','LineWidth',cfg.LineWidth)

plot([x3(3) x3(4)],[y3 y3],'k','LineWidth',cfg.LineWidth)

plot([x3(3) x3(3)],[y3 y4],'k','LineWidth',cfg.LineWidth)
plot([x3(4) x3(4)],[y3 y4],'k','LineWidth',cfg.LineWidth)

text(x1,y1,'*','FontSize',cfg.FontSize,'HorizontalAlignment','center')

% plot effect ex 2
subplot(1,5,5)
for iBar = 1:length(d4)
    h4(iBar) = bar(location(iBar),d4(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XLim',xlims,'YLim',ylims,'YTick',[])
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',cfg.FontSize,'LineWidth',1)
xlabel(xlab)
pbaspect(gca,[1 1 1])
box off; 

% plot the pval branches
plot([x1 x1],[y1 y2],'k','LineWidth',cfg.LineWidth)
plot([x2(1) x2(3)],[y2 y2],'k','LineWidth',cfg.LineWidth)
plot([x2(1) x2(1)],[y2 y3],'k','LineWidth',cfg.LineWidth)
plot([x2(3) x2(3)],[y2 y3],'k','LineWidth',cfg.LineWidth)
plot([x3(1) x3(2)],[y3 y3],'k','LineWidth',cfg.LineWidth)
plot([x3(3) x3(4)],[y3 y3],'k','LineWidth',cfg.LineWidth)
plot([x3(1) x3(1)],[y3 y4],'k','LineWidth',cfg.LineWidth)
plot([x3(2) x3(2)],[y3 y4],'k','LineWidth',cfg.LineWidth)
plot([x3(3) x3(3)],[y3 y4],'k','LineWidth',cfg.LineWidth)
plot([x3(4) x3(4)],[y3 y4],'k','LineWidth',cfg.LineWidth)

text(x1,y1,'*','FontSize',cfg.FontSize,'HorizontalAlignment','center')


%% save
cd(cfg.output_fd)
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-dpdf','-r300',cfg.output_fn);
print(gcf,'-depsc','-r300',cfg.output_fn);

cd(iWasHere)