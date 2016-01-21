%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%               Collect Arm Choice Behaviour and Plot                 %%%                
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Tmaze data, generates part of FIGURE_1 for paper

% ABOUT
% This script loads metadata and ExpKeys and collects behavioural data into
% a struct with the following organization:

%  behav.R042.food.L = mean nLeft trials on food day
%  behav.R042.food.R = mean nRight trials on food day
%  behav.R042.food.pLtrl = probably of going L on a food day across trials
%  behav.R042.water.L =  ...
%  . . .
%  behav.R042.water.pRtrl =  ^^      ^^      R on a water day    ^^   
%  behav.R044.food.L =  ...
%  . . .
%  and so on

% (food) and (water) refer to the restriction type, and (L) and (R) refer to the
% arms on the Tmaze.
% pLtrl fields, for example, contain the probability that the rat chose, L 
% (food reward) on a food restriction day for each trial.
% Note that "bad" and forced trials (as specified in the ExpKeys) are removed, 
% so not every 11th trial is the "real" 11th trial in terms of actual number of
% track traversals. These should probably be replaced with NaNs instead, but
% was done this way to preserve the original handling of the data in the cosyne 
% scripts. Some sessions have zero bad trials, and some sessions have two or 
% more. (A bad trial means the rat attempted to reverse directions and needed 
% to be interferred with by the researcher). Forced trials are quite
% common, and these are trials where the favoured arm was blocked and the
% rat had no choice but to go down the unfavoured arm.

% aacarey Sept 2015, based on thesis script which in turn was based on MvdM
%                    script MASTER_Behavior
% edit Jan 2015, turned two poster scripts into this script for the paper

clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO

cfg.writeDiary = 0; % keep a text file record of command window history
cfg.writeFiles = 1; % save data
cfg.output_fn = 'behavior'; % what to call the output

% Where to save the files
cfg.output_fd = [pwd,'\visuals'];  % file directory for saving data files

% Do you want to save the output?
cfg.writeFiles = 1; % 1 yes, 0 no

% What color scheme?
cfg.colormode = 'inventory4'; % see TmazeColors for colormode descriptions

% in the trial-by-choice subplot, do you want to plot the water inverted?
cfg.invertWater = 1; % 1 yes, 0 no
cfg.plotTrials = 1:15; % 1:15, plot trials 1 to 15 for the choice-by-trial figure

% do you want to display subpanel labels? (Like A, B, C)
cfg.showSubpanelLabels = 0; % 1 yes, 0 no

%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.plot);
    % and check has already been performed
else    
    % check that required things exist before continuing
    
    cfg.rats = {'R042','R044','R050','R064'};% recommended not to change the rats included, do all of them at the same time excepts when testing with the intention of deleting the test data
    cfg_temp = [];
    cfg_temp.requireExpKeys = 1;
    cfg_temp.ExpKeysFields = {'badTrials'};
    cfg_temp.requireMetadata = 1;
    cfg_temp.MetadataFields = {'taskvars'};
    if ~checkTmazeReqs(cfg_temp); return; end % make sure we have everything
    
end

%%

if cfg.writeDiary % save command window text
    warning off % I don't care, Navi
    cd(cfg.output_fd)
    diary(['behavior','.txt'])
    cd(cfg.output_fd)
    disp(date)
end

cfg.sessions = {'food','water'};
originalFolder = pwd;

cfg.rats = TmazeRats;
disp(['Rats: ',cfg.rats])

% collect data for all rats
% 'all' is the aggregator that combines data across rats
% 'behav' is the output struct, and contains indivual rat data and the
% combined rat data
all.nSessions = 0;
all.sessionChoice = []; % did they choose L or R
all.sessionType = []; % was it a food or a water day
all.nfood_chooseFood = 0; % the number of times a food restricted rat chose food reward
all.nfood_chooseWater = 0; % the number of times a food restricted rat chose water reward
all.nwater_chooseFood = 0;
all.nwater_chooseWater = 0;
behav.all.nLTrials = 0;
behav.all.nRTrials = 0;

for iRat = 1:length(cfg.rats)
    disp(' '); disp(cfg.rats(iRat))
    % give me a list of directories where I can find session data:
    cfg_temp = [];
    cfg_temp.rats = cfg.rats(iRat);
    fd = getTmazeDataPath(cfg_temp);
    
    %this holds the arm choices for each session separately
    sessionChoice = nan(length(fd),25); % length(fd) is same as number of sessions; 25 is more than the max trials done during any sessions
    sessionType = (1:length(fd))'; % will fill in session restriction type below: 1 = food, 2 = water
    all.nSessions = all.nSessions + length(fd);
    type = []; arm = [];
    nTrials = [];
    % go through each session and do the thing
    for iFD = 1:length(fd)
        cd(fd{iFD})
        
        LoadExpKeys
        LoadMetadata;
        
        sequence = metadata.taskvars.sequence;
        remove = sort(unique([ExpKeys.forcedTrials ExpKeys.badTrials]));
        sequence(remove) = [];
        nTrials = [nTrials; length(sequence)];
        
        restrictiontype = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
        type = [type;restrictiontype*ones(sum(strcmp(sequence,'L')==1)+sum(strcmp(sequence,'R')==1),1)];
        sessionType(iFD) = restrictiontype;
        
        arm = cat(1,arm,cat(1,ones(sum(strcmp(sequence,'L')==1),1)),2*ones(sum(strcmp(sequence,'R')==1),1));
        sequenceID = 1:length(sequence);
        for iTrial = 1:length(sequence)
            if restrictiontype == 1
                if strcmp(sequence(iTrial),'L')
                    sequenceID(iTrial) = 1; % for food on food day
                else
                    sequenceID(iTrial) = 0; % for water on food day
                end
            else
                if strcmp(sequence(iTrial),'L')
                    sequenceID(iTrial) = 0; % for food on water day
                else
                    sequenceID(iTrial) = 1; % for water on water day
                end
            end
        end
        
        sessionChoice(iFD,1:length(sequenceID)) = sequenceID;
    end
    
    %count food days and food day trials
    nfood_chooseFood = sum((arm == 1 & type == 1)); % arm 1 = left & type 1 = food
    nfood_chooseWater = sum((arm == 2 & type == 1));
    nFoodTrials = nfood_chooseFood + nfood_chooseWater;
    disp(['nfood ',num2str(nFoodTrials)])
    
    all.nfood_chooseFood = all.nfood_chooseFood + nfood_chooseFood;
    all.nfood_chooseWater = all.nfood_chooseWater + nfood_chooseWater;
    
    % ratio food trials
    food_chooseFoodR = nfood_chooseFood./nFoodTrials;
    food_chooseWaterR = nfood_chooseWater./nFoodTrials;
    
    disp(['Total food trials: ',num2str(nFoodTrials)])
    disp(['Chose food on food day: ',num2str(nfood_chooseFood)])
    disp(['Ratio: ',num2str(food_chooseFoodR)])
    disp(['Chose water on food day: ',num2str(nfood_chooseWater)])
    disp(['Ratio: ',num2str(food_chooseWaterR)])
    
    % count water days and water day trials
    nwater_chooseFood = sum((type == 2 & arm == 1));
    nwater_chooseWater = sum((type == 2 & arm == 2));
    nWaterTrials = nwater_chooseFood + nwater_chooseWater;
    
    disp(['nwater ',num2str(nWaterTrials)])
    all.nwater_chooseFood = all.nwater_chooseFood + nwater_chooseFood;
    all.nwater_chooseWater = all.nwater_chooseWater + nwater_chooseWater;
    
    % ratio water trials
    water_chooseFoodR = nwater_chooseFood./nWaterTrials;
    water_chooseWaterR = nwater_chooseWater./nWaterTrials;
    
    disp(' ')
    disp(['Total water trials: ',num2str(nWaterTrials)])
    disp(['Chose food on water day: ',num2str(nwater_chooseFood)])
    disp(['Ratio: ',num2str(water_chooseFoodR)])
    disp(['Chose water on water day: ',num2str(nwater_chooseWater)])
    disp(['Ratio: ',num2str(water_chooseWaterR)])
    
    disp(['Mean probability of choosing food on food days across trials: ',num2str(nanmean(sessionChoice(sessionType == 1,:)))])
    disp(['Mean probability of choosing water on water days across trials: ',num2str(nanmean(sessionChoice(sessionType == 2,:)))])
    
    % collect data into a struct
    
    % food restr data
    behav.(cfg.rats{iRat}).food.Lprob = food_chooseFoodR;
    behav.(cfg.rats{iRat}).food.Rprob = food_chooseWaterR;
    behav.(cfg.rats{iRat}).pLtrl = nanmean(sessionChoice(sessionType == 1,:));
    behav.(cfg.rats{iRat}).nLTrials = nFoodTrials;
    behav.all.nLTrials = behav.all.nLTrials + nFoodTrials;
    
    % water restr data
    behav.(cfg.rats{iRat}).water.Lprob = water_chooseFoodR;
    behav.(cfg.rats{iRat}).water.Rprob = water_chooseWaterR;
    behav.(cfg.rats{iRat}).pRtrl = nanmean(sessionChoice(sessionType == 2,:));
    behav.(cfg.rats{iRat}).nRTrials = nWaterTrials;
    behav.all.nRTrials = behav.all.nRTrials + nWaterTrials;
    
    % omg my script is too convoluted. collect these things for the
    % out-of-loop aggregation of all rat data
    all.sessionChoice = [all.sessionChoice; sessionChoice];
    all.sessionType = [all.sessionType; sessionType];
    
end

% aggregate data across rats and place in behav struct
behav.all.nSessions = all.nSessions;
behav.all.food.Lprob = all.nfood_chooseFood/behav.all.nLTrials;
behav.all.food.Rprob = all.nfood_chooseWater/behav.all.nLTrials;
behav.all.water.Lprob = all.nwater_chooseFood/behav.all.nRTrials;
behav.all.water.Rprob = all.nwater_chooseWater/behav.all.nRTrials;
behav.all.pLtrl = nanmean(all.sessionChoice(all.sessionType == 1,:));
behav.all.pRtrl = nanmean(all.sessionChoice(all.sessionType == 2,:));


if cfg.writeFiles
    cd(cfg.output_fd)
    save('behavior','behav')
end

if cfg.writeDiary; warning on; diary off; end

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

%~~~~~~~~~~~~ PLOT AVG BEHAVIOUR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

%~~~~~~~~~~~~ PLOT AVG TRIAL CHOICE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    plot(pRtrl,'Color',colors.all.w,'LineStyle','.','MarkerSize',20);
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


%~~~~~~~~~~~~ PLOT INDIVIDUAL BEHAVIOUR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
hold off

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

%% Finish up

if cfg.writeFiles
    cd(cfg.output_fd)
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-depsc','-r300',cfg.output_fn);
    print(gcf,'-dpdf','-r300',cfg.output_fn);
    %set(gcf,'Color','w')
    %img = getframe(gcf);
    %imwrite(img.cdata,[cfg.output_fn,'.png']);
end

warning on

cd(originalFolder)

disp('~~~~~~~~~~~ END OF BEHAVIOR RUN ~~~~~~~~~~~')