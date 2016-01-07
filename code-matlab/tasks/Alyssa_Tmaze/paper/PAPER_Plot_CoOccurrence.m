%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%          Plot Prerecord, Task, Postrecord CoOccurrence Data         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This takes in the coocP output from PAPER_CoOccurrence and generates a 
% multiplot of all co-occurrence data for either p0 or p4 depending on
% which one is selected in cfg.whichP.
% At the time this was written, the script assumes you have the following
% folder organization:
% TmazePaper\data
% TmazePaper\visuals
%
% Make sure you have your current directory set to TmazePaper (or whatever 
% the main folder is called), and the script should run as expected 
% supposing that the following data has been generated:
%
% coOccurrence_prerecord.mat
% coOccurrence_task.mat *OR* coOccurrence_allITI.mat
% coOccurence_postrecord.mat
%
% The above data is generated from different passes through the script
% PAPER_CoOccurrence with the only difference being cfg.whichEvents.
%
% aacarey Nov 2015, Dec 2015, from POSTER version of this script

clear % note that the base workspace is cleared

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% what to plot?
cfg.whichP = 'p4m'; % 'p0' or 'p4' (single cell or zscore coactivity), can also do p1, p2, p3, p5 now

cfg.whichTask = 'task'; % 'task' for all events occurring during the task 
% (both track and waiting platform SWRs) or 'allITI' for all events 
% occurring during the intertrial intervals when the rat was sitting on the
% waiting platform.

% some appearance things
cfg.colormode = 'inventory3';
FontSize = 8;

% where to save the figure?
cfg.input_fd = [pwd,'\data'];
cfg.output_fd = [pwd,'\visuals'];

% Controls the amount of whitespace above and below the y data min and max
% with respect to the y ticks and labels. This is an "at least" value
tickBuffer = 0.1; % at least this much buffer space, as a proportion of the data extremes


%%
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)

% get prerecord co-occur data
pre = load('coOccurrence_prerecord.mat','coocP');

% get task or allITI co-occur data
switch cfg.whichTask
    case 'task'
        task = load('coOccurrence_task.mat','coocP');
    case 'allITI'
        task = load('coOccurrence_allITI.mat','coocP');
end

% get postrecord co-occur data
post = load('coOccurrence_postrecord.mat','coocP');

% some things that vary depending on which p is chosen:
switch cfg.whichP
    case 'p4'
        ylab = {'SWR coactivation'; 'Z-score'};
        tickDiff_all = 0.25;
        tickDiff_single = 1;
        cfg.output_fn = 'p4Multi';
        
    case 'p4m'
        ylab = {'dlfkgjlkjdh'; 'Z-score'};
        tickDiff_all = 5*10^-17;
        tickDiff_single = 5*10^-17;
        cfg.output_fn = 'p4mMulti';
        
    case 'p5'
        ylab = {'Randomly shuffled'; 'stuff'};
        tickDiff_all = 0.05;
        tickDiff_single = 0.1;
        cfg.output_fn = 'p5Multi';
        
    case 'p3'
        ylab = {'Cell pair'; 'joint probability'};
        tickDiff_all = 0.01;
        tickDiff_single = 0.1;
        cfg.output_fn = 'p3Multi';
        
    case 'p2'
        ylab = {'Cell pair'; 'conditional probability'};
        tickDiff_all = 0.05;
        tickDiff_single = 0.1;
        cfg.output_fn = 'p2Multi';
        
    case 'p1'
        ylab = {'Cell pair'; 'whatever this is'};
        tickDiff_all = 0.01;
        tickDiff_single = 0.1;
        cfg.output_fn = 'p1Multi';
        
    case 'p0'
        ylab = {'Proportion of';'SWRs active'};
        tickDiff_all = 0.05;
        tickDiff_single = 0.1;
        cfg.output_fn = 'p0Multi';
        
    otherwise
        error('Unrecognized cfg.whichP')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                Section for combined data plots                      %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here's where we plot the overall data (across 4 rats)

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    

figure; hold on

%~~~~~~~~~~~~~~~~ PRERECORD SUBPLOT (ALL RATS) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(4,6,[1 2 7 8])
d = [nanmean(pre.coocP.(rats{iRat}).foodL.(cfg.whichP)) nanmean(pre.coocP.(rats{iRat}).foodR.(cfg.whichP)) nanmean(pre.coocP.(rats{iRat}).waterL.(cfg.whichP)) nanmean(pre.coocP.(rats{iRat}).waterR.(cfg.whichP))];
for iBar = 1:length(d)
    bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

h_all(1) = gca;
mindata_all(1) = min(d); % for y limits assignment, get lowest value
maxdata_all(1) = max(d); % for y limits, get highest value

set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YTick',[])
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

%~~~~~~~~~~~~~~~~~~~ TASK SUBPLOT (ALL RATS) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(4,6,[3 4 9 10])
d = [nanmean(task.coocP.(rats{iRat}).foodL.(cfg.whichP)) nanmean(task.coocP.(rats{iRat}).foodR.(cfg.whichP)) nanmean(task.coocP.(rats{iRat}).waterL.(cfg.whichP)) nanmean(task.coocP.(rats{iRat}).waterR.(cfg.whichP))];
for iBar = 1:length(d)
    bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

h_all(2) = gca; % get axes handle
mindata_all(2) = min(d); % for y limits assignment, get lowest value
maxdata_all(2) = max(d); % for y limits, get highest value

set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('TASK')
box off
set(gca,'Layer','top')
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)

%~~~~~~~~~~~~~~~~ POSTRECORD SUBPLOT (ALL RATS) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
subplot(4,6,[5 6 11 12]) % for all rats POST
d = [nanmean(post.coocP.(rats{iRat}).foodL.(cfg.whichP)) nanmean(post.coocP.(rats{iRat}).foodR.(cfg.whichP)) nanmean(post.coocP.(rats{iRat}).waterL.(cfg.whichP)) nanmean(post.coocP.(rats{iRat}).waterR.(cfg.whichP))];
for iBar = 1:length(d)
    bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

h_all(3) = gca;
mindata_all(3) = min(d); % for y limits assignment, get lowest value
maxdata_all(3) = max(d); % for y limits, get highest value

set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('POSTRECORD')
box off
set(gca,'Layer','top')
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)

%~~~~~~~~~ some auto formatting of y axis stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% roughly controls how many ticks labels there are...not worth time to make it perfect
nTickLabels_all = 6;

% set y ticks and limits for combined data subplots (the big ones at the top)
mindata_all = min(mindata_all); mindata_all = min([0 floor((mindata_all + mindata_all*tickBuffer)/tickDiff_all)*tickDiff_all]);
maxdata_all = max(maxdata_all); maxdata_all = ceil((maxdata_all + maxdata_all*tickBuffer)/tickDiff_all)*tickDiff_all;

%stepSize = round((diff([maxdata_all mindata_all])/5)/5)*5; % wut

yticks_all = mindata_all:tickDiff_all:maxdata_all;
yTickStep = round(length(yticks_all)/(nTickLabels_all));
yticks_all = yticks_all(1:yTickStep:length(yticks_all));
%yticks_all(end) = yticks_all(end) + diff([yticks_all(end-1) yticks_all(end)]);

set(h_all(1),'YTick',yticks_all)
set(h_all,'YLim',[yticks_all(1),yticks_all(end)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                 Section for single rat plots                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here's where we plot the data for each rat individually

rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.65; yName = 0.85; % controls where to place rat name text

axisNumber = 1; % tracks which subplot we're working with: main purpose is to make auto assignment of y limits and y ticks possible (see code in loop below)
% preallocate space for tracking axes handles and data value extremes (just because I hate seeing orange stuff in my scripts :P)
h_single = nan(12,1); % "empty" container for axes handles
mindata_single = h_single; % empty container for data value lower extremes
maxdata_single = h_single; % empty container for data value higher extremes

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4]; % tracks which spot the subplot will occupy
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    %~~~~~~~~~~~~ PRERECORD SUBPLOT (SINGLE RAT) ~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    d = [nanmean(pre.coocP.(rats{iRat}).foodL.(cfg.whichP)) nanmean(pre.coocP.(rats{iRat}).foodR.(cfg.whichP)) nanmean(pre.coocP.(rats{iRat}).waterL.(cfg.whichP)) nanmean(pre.coocP.(rats{iRat}).waterR.(cfg.whichP))];
    for iBar = 1:length(d)
        bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    if position(1) == 13 || position(1) == 19
        ylabel(ylab,'FontSize',FontSize)
    else
        yticks = [];
    end
    
    % track axes and min and max data values
    h_single(axisNumber) = gca;
    axisNumber = axisNumber+1;
    mindata_single(axisNumber) = min(d); % for y limits assignment, get lowest value
    maxdata_single(axisNumber) = max(d); % for y limits, get highest value
    
    % set some axes properties
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    %~~~~~~~~~~~~ TASK SUBPLOT (SINGLE RAT) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subplot(4,6,position(2))
    d = [nanmean(task.coocP.(rats{iRat}).foodL.(cfg.whichP)) nanmean(task.coocP.(rats{iRat}).foodR.(cfg.whichP)) nanmean(task.coocP.(rats{iRat}).waterL.(cfg.whichP)) nanmean(task.coocP.(rats{iRat}).waterR.(cfg.whichP))];
    for iBar = 1:length(d)
        bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    % track axes and min and data extremes
    h_single(axisNumber) = gca;
    axisNumber = axisNumber+1;
    mindata_single(axisNumber) = min(d); % for y limits assignment, get lowest value
    maxdata_single(axisNumber) = max(d); % for y limits, get highest value
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    %~~~~~~~~~~~~ POSTRECORD SUBPLOT (SINGLE RAT) ~~~~~~~~~~~~~~~~~~~~~~~~~
    subplot(4,6,position(3))
    d = [nanmean(post.coocP.(rats{iRat}).foodL.(cfg.whichP)) nanmean(post.coocP.(rats{iRat}).foodR.(cfg.whichP)) nanmean(post.coocP.(rats{iRat}).waterL.(cfg.whichP)) nanmean(post.coocP.(rats{iRat}).waterR.(cfg.whichP))];
    for iBar = 1:length(d)
        bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    % track axes and data extremes
    h_single(axisNumber) = gca;
    axisNumber = axisNumber+1;
    mindata_single(axisNumber) = min(d); % for y limits assignment, get lowest value
    maxdata_single(axisNumber) = max(d); % for y limits, get highest value
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

%~~~~~~~~~ some auto formatting of y axis stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% kind of controls how many ticks labels there are...not worth time to make it perfect
nTickLabels_single = 4;

% set y ticks and limits for combined data subplots (the big ones at the top)
mindata_single = min(mindata_single); mindata_single = min([0 round((mindata_single + mindata_single*tickBuffer)/tickDiff_single)*tickDiff_single]);
maxdata_single = max(maxdata_single); maxdata_single = round((maxdata_single + maxdata_single*tickBuffer)/tickDiff_single)*tickDiff_single;

yticks_single = mindata_single:tickDiff_single:maxdata_single;
yTickStep = round(length(yticks_single)/(nTickLabels_single));
yticks_single = yticks_single(1:yTickStep:length(yticks_single));
%yticks_single(end) = yticks_single(end) + diff([yticks_single(end-1) yticks_single(end)]);

set(h_single(1),'YTick',yticks_single)
set(h_single(7),'YTick',yticks_single)
set(h_single,'YLim',[yticks_single(1),yticks_single(end)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        SAVE SECTION                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(cfg.output_fd); 
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-depsc','-r300',cfg.output_fn);
print(gcf,'-dpdf','-r300',cfg.output_fn);
cd(originalFolder)