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

clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% what to plot?
cfg.whichP = 'p4'; % 'p0' or 'p4' (single cell or zscore coactivity), can also do p1, p2, p3, p5 now

cfg.whichTask = 'task'; % 'task' for all events occurring during the task 
% (both track and waiting platform SWRs) or 'allITI' for all events 
% occurring during the intertrial intervals when the rat was sitting on the
% waiting platform.

% some appearance things
cfg.colormode = 'inventory4';
FontSize = 8;

% where to save the figure?
cfg.input_fd = [pwd,'\data'];
cfg.output_fd = [pwd,'\visuals'];

% Controls the amount of whitespace above and below the y data min and max
% with respect to the y ticks and labels. This is an "at least" value
tickBuffer = 0.1; % at least this much buffer space, as a proportion of the data extremes

%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.plot);
end

%%
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)

% get prerecord co-occur data
data.prerecord = loadpop('coOccurrence_prerecord.mat','coocP');

% get equal behaviorITI data
data.equalBehaviorITI = loadpop('coOccurrence_equalBehaviorITI.mat','coocP');

% get task or allITI co-occur data
switch cfg.whichTask
    case 'task'
        data.task = loadpop('coOccurrence_task.mat','coocP');
    case 'allITI'
        data.allITI = loadpop('coOccurrence_allITI.mat','coocP');
end

% get postrecord co-occur data
data.postrecord = loadpop('coOccurrence_postrecord.mat','coocP');

% some things that vary depending on which p is chosen:
switch cfg.whichP
    case 'p4'
        ylab = {'SWR coactivation'; 'Z-score'};
        tickDiff_all = 0.25;
        tickDiff_single = 0.1;
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
 
rats = ['all',TmazeRats];

epochs = fieldnames(data);

position_all = [1 2 9 10; 3 4 11 12; 5 6 13 14; 7 8 15 16];
%start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
start = [17 18 25 26];
jhgjhg = [0 2 4 6];


axisNumber = 1; % tracks which subplot we're working with: main purpose is to make auto assignment of y limits and y ticks possible (see code in loop below)
% preallocate space for tracking axes handles and data value extremes (just because I hate seeing orange stuff in my scripts :P)
h_single = nan(12,1); % "empty" container for axes handles
mindata_single = h_single; % empty container for data value lower extremes
maxdata_single = h_single; % empty container for data value higher extremes


figure; hold on

for iEpoch = 1:length(epochs)
    for iRat = 1:length(rats)
        switch rats{iRat}
            case 'all'
                position = position_all(iEpoch,:);
                txt = 'all rats';
            otherwise
                position = start(iRat-1) + jhgjhg(iEpoch); % tracks which spot the subplot will occupy
                txt = rats{iRat};
        end
        
        col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
        
        subplot(4,8,position)
        d = [nanmean(data.(epochs{iEpoch}).(rats{iRat}).foodL.(cfg.whichP)) nanmean(data.(epochs{iEpoch}).(rats{iRat}).foodR.(cfg.whichP)) nanmean(data.(epochs{iEpoch}).(rats{iRat}).waterL.(cfg.whichP)) nanmean(data.(epochs{iEpoch}).(rats{iRat}).waterR.(cfg.whichP))];
        for iBar = 1:length(d)
            bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
            hold on
        end
        
        if iRat == 1
            h_all(iEpoch) = gca;
            mindata_all(iEpoch) = min(d); % for y limits assignment, get lowest value
            maxdata_all(iEpoch) = max(d); % for y limits, get highest value
            
            % add title
            switch epochs{iEpoch}
                case 'prerecord'
                    titl = 'PRERECORD';
                case 'allITI'
                    titl = 'INTERTRIAL INTERVALS';
                case 'task'
                    titl = 'TASK';
                case 'equalBehaviorITI'
                    titl = 'EQUAL BEHAVIOR ITI';
                case 'postrecord'
                    titl = 'POSTRECORD';
                otherwise
                    titl= 'unassigned epoch name!';
            end
            title(titl,'FontSize',FontSize)
            
            xlabel('  food                  water','FontSize',FontSize)
            
            text(0.58,0.9,'all rats','Units','normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','right',...
                'FontSize',FontSize)
            
        else
            % track axes and min and max data values
            h_single(axisNumber) = gca;
            axisNumber = axisNumber+1;
            mindata_single(axisNumber) = min(d); % for y limits assignment, get lowest value
            maxdata_single(axisNumber) = max(d); % for y limits, get highest value
            
            xlabel('  food   water','FontSize',FontSize)
            
            text(0.65,0.85,rats{iRat},'Units','normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','right',...
                'FontSize',FontSize)
        end
        
        if isequal(position,[1 2 9 10]) || position(1) == 17 || position(1) == 25
            ylabel(ylab,'FontSize',FontSize)
        end
        
        box off
        set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
        %set(gca,'YTick',[])
        set(gca,'Layer','top')
        
    end % of rats
end % of epochs

%~~~~~~~~~ some auto formatting of y axis stuff for avg rat plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

if ~ isempty(yticks_all)
    set(h_all,'YTick',[])
    set(h_all(1),'YTick',yticks_all)
    set(h_all,'YLim',[yticks_all(1),yticks_all(end)])
end
   
%~~~~~~~~~ some auto formatting of y axis stuff for single rat plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~

% kind of controls how many ticks labels there are...not worth time to make it perfect
nTickLabels_single = 4;

% set y ticks and limits for combined data subplots (the big ones at the top)
mindata_single = min(mindata_single); 
mindata_single = min([0 round((mindata_single + mindata_single*tickBuffer)/tickDiff_single)*tickDiff_single]);
maxdata_single = max(maxdata_single); 
maxdata_single = round((maxdata_single + maxdata_single*tickBuffer)/tickDiff_single)*tickDiff_single;

yticks_single = mindata_single:tickDiff_single:maxdata_single;
yTickStep = round(length(yticks_single)/(nTickLabels_single));
yticks_single = yticks_single(1:yTickStep:length(yticks_single));
yticks_single = [yticks_single , yticks_single(end) + diff([yticks_single(end-1) yticks_single(end)])];

set(h_single(1),'YTick',yticks_single)
set(h_single(3),'YTick',yticks_single)
set(h_single,'YLim',[yticks_single(1),yticks_single(end)])

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        SAVE SECTION                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(cfg.output_fd); 
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-depsc','-r300',cfg.output_fn);
print(gcf,'-dpdf','-r300',cfg.output_fn);
cd(originalFolder)

disp('~~~~~~~ END OF PLOT RUN ~~~~~~~')