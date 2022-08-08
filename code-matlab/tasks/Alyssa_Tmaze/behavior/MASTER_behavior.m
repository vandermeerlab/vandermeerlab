%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASTER_behavior.m %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generates Figure 1b in Carey & van der Meer (2015) and associated stats
%
% NOTE: getTmazeDataPath() will need to be configured correctly for your
% machine. This is the function that needs to know where your data files
% are.

%% get sessions we can use
cfg = []; cfg.rats = {'R050'}; cfg.writeFig = 0;
fd = getTmazeDataPath(cfg);

nSessions = length(fd);
nMaxTrials = 25;

% init vars
out = [];
out.session_id = cell(nSessions,1);
out.rat_id = cell(nSessions,1);
out.type = cell(nSessions,1);

out.nTrials = nan(nSessions,1);
out.nWater = nan(nSessions,1);
out.nFood = nan(nSessions,1);

out.fullChoice = nan(nSessions,nMaxTrials);

curr_pwd = pwd; % store curr dir so can revert to later

%% collect data
for iSession = 1:nSessions
    
    cd(fd{iSession});
    out.session_id{iSession} = pwd; % should only store SessionID
        
    %
    load(FindFile('*metadata.mat'));
    run(FindFile('*keys.m'));
    
    %
    out.type{iSession} = ExpKeys.RestrictionType;
    
    %
    this_sequence = metadata.taskvars.sequence;
    %this_sequence(ExpKeys.badTrials) = []; %@aacarey, what is best to do with these?
    this_sequence(ExpKeys.forcedTrials) = []; % exclude forced trials
    
    %
    out.nTrials(iSession) = length(this_sequence);
    
    
    switch ExpKeys.Layout
        case 'foodLeft'
            foodLoc = 'L'; waterLoc = 'R';
        case 'foodRight'
            foodLoc = 'R'; waterLoc = 'L';
    end
        
    this_sequence = cellfun(@(x) strcmp(x,foodLoc),this_sequence); % now food is 1, water 0
    
    out.nWater(iSession) = sum(this_sequence == 0);
    out.nFood(iSession) = sum(this_sequence == 1);    
    
    out.fullChoice(iSession,1:length(this_sequence)) = this_sequence;
    
end % over sessions

%% process data into format useful for plotting and stats
type = cellfun(@(x) strcmp(x,'food'),out.type); % food is 1, water is 0

% counts
nFoodTrials = sum(out.nTrials(type == 1));
food_chooseFoodN = sum(out.nFood(type == 1));
food_chooseWaterN = sum(out.nWater(type == 1));

% ratio as total trials
food_chooseFoodR = food_chooseFoodN./nFoodTrials;
food_chooseWaterR = food_chooseWaterN./nFoodTrials;

% counts
nWaterTrials = sum(out.nTrials(type == 0));
water_chooseFoodN = sum(out.nFood(type == 0));
water_chooseWaterN = sum(out.nWater(type == 0));

% ratio as total trials
water_chooseFoodR = water_chooseFoodN./nWaterTrials;
water_chooseWaterR = water_chooseWaterN./nWaterTrials;

%% make Figure 1b
figure(1);

d = [food_chooseFoodR food_chooseWaterR water_chooseFoodR water_chooseWaterR]; % from excel sheet
xi = [1 2 4 5];
cl = 'rrbb';

for iB = 1:length(d)
   
    h(iB) = bar(xi(iB),d(iB)); set(h(iB),'FaceColor',cl(iB),'EdgeColor','none');
    hold on;
    
end

set(gca,'XTick',xi,'XTickLabel',{'food','water','food','water'}, ...
    'YLim',[0 1],'LineWidth',1,'YTick',0:0.25:1,'FontSize',18);

ylabel('p(choice)');
xlabel('(food restricted)               (water restricted)');
title(sprintf('%d rats, %d sessions, %d food trials, %d water trials',length(cfg.rats),nSessions,nFoodTrials,nWaterTrials));

box off;

%% generate stats for Figure 1b

% to be done

%% make Figure S1a
figure(2); clear h;

% for type food is 1, water is 0
h(1) = plot(nanmean(out.fullChoice(type == 0,:)),'b'); % water restrict
hold on; 
plot(nanmean(out.fullChoice(type == 0,:)),'.b','MarkerSize',20);

h(2) = plot(nanmean(out.fullChoice(type == 1,:)),'r'); % food restrict
plot(nanmean(out.fullChoice(type == 1,:)),'.r','MarkerSize',20);

set(gca,'XTick',0:5:25, ...
    'YLim',[0 1.01],'LineWidth',1,'YTick',0:0.25:1,'FontSize',18);
box off;

xlabel('trial number'); ylabel('p(choose food)');

legend(h,{'water-restricted','food-restricted'},'Location','Southeast');
legend boxoff;

%% write
cd(curr_pwd);
if cfg.writeFig
   
    base_figname = 'MASTER_behavior_fig';
    
    for iF = 1:2
       
        figure(iF);
        print(gcf,'-dpng','-r300',cat(2,base_figname,num2str(iF),'.png'));
        print(gcf,'-dill',cat(2,base_figname,num2str(iF),'.ai'));
        
    end
    
    
end
