%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                         Plot Tmaze Behavior                         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the script 'Behavior_GenData' before running this script. The scripts 
% exist separately to simplify / speed up the process of plotting data and
% making changes to the appearance of figures. 

% Generates a multipanel figure containing behavior for individual rats and
% all rats together. At the time this was written, the output figure was
% generated so that no postprocessing steps were necessary.
% The behavioral data includes free choice trials only (no forced/blocked
% trials), and bad trials (eg runbacks) have been removed.

% The settings are such that the png / pdf files are properly arranged. Do
% not consider the figure MATLAB shows to be the final appearance. Instead,
% look at the saved file.

% This script is modified from aacarey's thesis Behavior script, which
% was modified from MvdM's MASTER_behavior.

% aacarey Sept 2015

clear

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% no cfg.rats field; must plot all rats with this script

% do you want to save the output?
cfg.writeFiles = 1; % 1 yes, 0 no

cfg.colormode = 'inventory3'; % see TmazeColors for colormode descriptions
% colormodes:
%    'inventory1': indiv rats are inventory colors; combined is grey
%    'inventory2': indiv rats are inventory; combined is red and blue
%    'rb': everything in red and blue
%    'grey': everything in grey

% what to load and where to find it
cfg.input_fn = 'behavior'; % the script knows the proper file extension (.mat) so don't put it here
cfg.input_fd = 'E:\Documents\TmazePaper\data';

% what to call the output figure and where to save it
cfg.output_fn = 'behaviorMulti'; % it saves as a .png
cfg.output_fd = 'E:\Documents\TmazePaper\visuals';

% in the trial-by-choice subplot, do you want to plot the water inverted?
cfg.invertWater = 1; % 1 yes, 0 no
cfg.plotTrials = 1:15; % 1:15, plot trials 1 to 15 for the choice-by-trial figure

% do you want to display subpanel labels? (Like A, B, C)
cfg.showSubpanelLabels = 0; % 1 yes, 0 no


%% Load the data

iWasHere = pwd; % remember where we started

cd(cfg.input_fd)
load(FindFile([cfg.input_fn,'.mat'])); % loads struct called 'behav'

% get colors to plot with
colors = TmazeColors(cfg.colormode); % colors is struct with food day and water day 
% colors for each rat. In some cases they might all be the same, but others 
% they are all different. this allows flexible plotting of different colors, if wanted


warning off

%% make plot

rats = {'R042','R044','R050','R064'}; % do not change {'R042','R044','R050','R064'}

location = [1 2 4 5]; % where to place the bar

% get the accumulated data subplotted first 

FontSize = 10; 

figure; hold on

subplot(3,4,[1 2 5 6])
    % what's it say down there? d = ['probability of going left on food day' 'probability of going right on food day' 'probability of going left on water day' 'probability of going right on water day'];
    d = [behav.all.food.Lprob behav.all.food.Rprob behav.all.water.Lprob behav.all.water.Rprob];
    % the above evaluates to (ex:)  behav.R042.food.Lprob ... which is the probability the rat went left on a food restriction day
    col = {colors.all.f colors.all.f colors.all.w colors.all.w}; % {foodColor foodColor waterColor waterColor}
    
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        %errorbar(location(iBar),all.meanLatency(iBar),all.ysem(iBar),'k','LineStyle','none');
        hold on
    end
    
    ylabel('p(choice)','FontSize',FontSize); set(gca,'XLim',[0 location(4)+1])
    set(gca,'YLim',[0 1],'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'YLim',[0 1],'LineWidth',1,'YTick',0:0.25:1);
    xlabel('       food restricted               water restricted      ','FontSize',FontSize)
    
    title(['all rats: ',num2str(behav.all.nLTrials),' left trials, ',num2str(behav.all.nRTrials),' right trials']);
    box off
    set(gca,'Layer','top')

    % this controls the position of the subpanel label for "all rats"
    if cfg.showSubpanelLabels
    txt = 'A';
    text(-.1,0.95,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',24)
    end
    
    % this controls the position of the "all rats" label; number values
    % need to be adjusted if font size is changed
    %txt = 'all rats';
    %text(0.56,0.9,txt,'Units','normalized',...
    %    'VerticalAlignment','bottom',...
    %    'HorizontalAlignment','right',...
    %    'FontSize',FontSize)
    
subplot(3,4,[3 4 7 8])
    
    pLtrl = behav.all.pLtrl(cfg.plotTrials);
    pRtrl = behav.all.pRtrl(cfg.plotTrials);
    
    if cfg.invertWater % control some plotting visuals (distinctiveness)
        foodstyle = '.';
        foodsize = 20;
    else
        foodstyle = 'x';
        foodsize = 10;
        
    end
    
    % plot food choice on food days
    h(1) = plot(pLtrl,'Color',colors.all.f,'LineWidth',2); 
    hold on;
    plot(pLtrl,'Color',colors.all.f,'LineStyle',foodstyle,'MarkerSize',foodsize);
    
    % plot water choice on water days
    if cfg.invertWater
        ugh = [pRtrl; ones(size(pRtrl))];
        pRtrl = abs(diff(ugh,1));
        h(2) = plot(pRtrl,'Color',colors.all.w,'LineWidth',2);
        %plot(pRtrl,'Color',colors.all.w,'LineStyle','.','MarkerSize',20);
        ylims = [0 1];
        ylabel('p(choose food)')
    else
        h(2) = plot(pRtrl,'Color',colors.all.w,'LineWidth',2);
        plot(pRtrl,'Color',colors.all.w,'LineStyle','.','MarkerSize',20);
        ylims = [0.25 1];
    end
    
    set(gca,'XTick',0:5:25, ...
        'YLim',ylims,'XLim',[0 length(cfg.plotTrials)+2],'LineWidth',1,'YTick',0.25:0.25:1,'FontSize',FontSize);
    box off;
    set(gca,'Layer','top')
    
    %title(['All rats, 24 sessions, ',num2str(behav.all.nFood),' left trials, ',num2str(behav.all.nFood), ' right trials']);
    xlabel('trial number','FontSize',FontSize); %ylabel('p(choice)');
    
    legend(h,{'food restricted','water restricted'},'Location','East','FontSize',FontSize);
    legend boxoff;
    
    % this controls the position of the subpanel label for the
    % choice-by-trial plot
    if cfg.showSubpanelLabels
    txt = 'B';
    text(-.1,0.95,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',24)
    end
    

    % work on individual rat plots
    
    position = [9,10,11,12];
    
for iRat = 1:length(rats)
    
    % what's it say down there? d = ['probability of going left on food day' 'probability of going right on food day' 'probability of going left on water day' 'probability of going right on water day'];
    d = [behav.(rats{iRat}).food.Lprob behav.(rats{iRat}).food.Rprob behav.(rats{iRat}).water.Lprob behav.(rats{iRat}).water.Rprob];
    % the above evaluates to (ex:)  behav.R042.food.Lprob ... which is the probability the rat went left on a food restriction day
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    subplot(3,4,position(iRat))
    
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    % control y label...show it for leftmost plot only
    if position(iRat) == 9
        ylabel('p(choice)','FontSize',FontSize);
        
        % this controls the position of the subpanel label for the
        % individual rat plots
        if cfg.showSubpanelLabels
        txt = 'C';
        text(-.24,0.95,txt,'Units','normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','right',...
            'FontSize',24)
        end
        
    end
    
    set(gca,'XLim',[0 location(4)+1])
    set(gca,'YLim',[0 1],'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'YLim',[0 1],'LineWidth',1,'YTick',0:0.25:1);
    xlabel('  food          water ', 'FontSize',FontSize)
    
    %title([rats{iRat},', 6 sessions, ',num2str(behav.(rats{iRat}).nFood),' left trials, ',num2str(behav.(rats{iRat}).nWater), ' right trials']);
    box off
    set(gca,'Layer','top')
    
    % add rat id inside of figure axes; number values
    % need to be adjusted if font size is changed
    txt = rats{iRat};
    text(0.63,0.8,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

%maximize

if cfg.writeFiles
    cd(cfg.output_fd)
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']); % this doesn't save with the dimension i see on screen
    print(gcf,'-depsc','-r300',cfg.output_fn);
    print(gcf,'-dpdf','-r300',cfg.output_fn);
    %set(gcf,'Color','w')
    %img = getframe(gcf);
    %imwrite(img.cdata,[cfg.output_fn,'.png']);
end

cd(iWasHere)
warning on