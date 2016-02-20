% batch script to collect and plot co-occurrence data

% collect co-occurrence data so it can be easily put into a multiplot 
% SfN poster *rushcode*

%%

cfg.output_fd = 'E:\Documents\TmazePaper\data';
cfg.writeDiary = 1; % keep a record of command window history in a text file
cfg.printFig = 0; % save figure (note: figure still plots)

cfg.whichEvents = 'all'; % 'all','prerecord','postrecord','taskrest','task'

cfg.saveData = 0;


cfg.colormode = 'inventory2';

%% do not change these configs.

cfg.what = 'SWRz'; % 'SWRz', use only the data for co-occurrence during candidate events
cfg.rats = {'R042','R044','R050','R064'};
switch cfg.whichEvents
    case 'all'
        cfg.prefix = 'CC_ALL_';
    case 'prerecord'
        cfg.prefix = 'CC_PRE_';
    case 'postrecord'
        cfg.prefix = 'CC_POST_';
        
    case 'taskrest'
        cfg.prefix = 'CC_TASKREST_';
        
    case 'task'
        cfg.prefix = 'CC_TASK_';
    otherwise
        error('Unrecognized cfg.whichEvents')
end

cfg.requireCandidates = 1;
fd = getTmazeDataPath(cfg);

cfg.sessions = {'food','water'};
cfg.arms = {'left','right'}; % to make ALL_obs_p0.arm


originalFolder = pwd;
if cfg.writeDiary % save command window text
    cd(cfg.output_fd)
    diary([cfg.prefix,'stats_rush.txt'])
    cd(originalFolder)
    disp(date)
end

disp(['Script name: ',mfilename])
disp(' ')


switch cfg.whichEvents
    case 'all'
        disp('** ALL EVENTS **')
        
    case 'prerecord'
        disp('** PRERECORD EVENTS **')
        
    case 'postrecord'
        disp('** POSTRECORD EVENTS **')
        
    case 'taskrest'
        disp('** TASKREST EVENTS **')
        
    case 'task'
        disp('** TASK EVENTS **')
end
disp(' ')
disp('You have selected: ')
disp(cfg.rats)
disp(' ')

%% GET COMBINED DATA
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%                                                               %%%')
disp('%%%             Combined co-occurrence data (all rats):           %%%')
disp('%%%                                                               %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% init vars to place data into

ALL_obs_p0.data = []; % the actual observed activation probabilities
ALL_obs_p0.arm = []; % left (1), right (2)
ALL_obs_p0.type = []; % restriction type (food/water)
ALL_obs_p0.sess = []; % session ID in case we want to restrict later

ALL_obs_p4.data = []; % the actual observed activation probabilities
ALL_obs_p4.arm = []; % left (1), right (2)
ALL_obs_p4.type = []; % restriction type (food/water)
ALL_obs_p4.sess = []; % session ID in case we want to restrict later

% the data that is loaded in the next section contains

%out.run.full_cc = c_crossed;
%out.run.left = all_left;
%out.run.right = all_right;
%out.run.crossed = all_crossed;
%  data from in-field coactivity (unused)

%out.nLeft 
%  the number of left-only cells

%out.nRight 
%  the number of right-only cells

%out.TCleft = TC_left; out.TCright = TC_right;

%out.SWRz.full_cc = c_crossed;
%  all co-occurrence pvals made from a Q-matrix composed of both left and
%  right track cells

%out.SWRz.left = all_left;
%  p4 for left-only pairs

%out.SWRz.right = all_right;
%  p4 for right-only pairs

%out.SWRz.crossed = all_crossed;
%  p4 for left-right pairs

%%

for iFD = 1:length(fd)
    
   cd(fd{iFD}); 
   LoadExpKeys;
   
   cd([pwd,'\','files']);
   this_file = FindFile(cat(2,cfg.prefix,'*coOccur_data.mat'));
   load(this_file);
   
   %
   this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));

   % collect data
   
   % get data to process from loaded file
   this_data = out.(cfg.what); %eval(cat(2,'out.',cfg.what,';'));
   nL = size(this_data.left,1); % number of left cells
   nR = size(this_data.right,1); % number of right cells
   
   tempL0 = this_data.full_cc.p0(1:nL); tempR0 = this_data.full_cc.p0(nL+1:end); 
   if size(tempL0,2) == 0, tempL0 = tempL0'; end; if size(tempR0,2) == 0, tempR0 = tempR0'; end
   temp0 = cat(1,tempL0,tempR0);
   ALL_obs_p0.data = cat(1,ALL_obs_p0.data,temp0); % the actual observed activation probabilities
   ALL_obs_p0.arm = cat(1,ALL_obs_p0.arm,cat(1,ones(length(tempL0),1),2*ones(length(tempR0),1))); % this make an identity array
   % thing that contains 1's in places corresponding to left cell data (probabilities), and 2's for right cell data
   ALL_obs_p0.type = cat(1,ALL_obs_p0.type,this_session_type*ones(size(temp0)));
   ALL_obs_p0.sess = cat(1,ALL_obs_p0.sess,iFD*ones(size(temp0)));
   
   tempL4 = this_data.left(:); tempR4 = this_data.right(:); temp4 = cat(1,tempL4,tempR4);
   ALL_obs_p4.data = cat(1,ALL_obs_p4.data,temp4); % the actual observed activation probabilities
   ALL_obs_p4.arm = cat(1,ALL_obs_p4.arm,cat(1,ones(length(tempL4),1),2*ones(length(tempR4),1)));
   ALL_obs_p4.type = cat(1,ALL_obs_p4.type,this_session_type*ones(size(temp4)));
   ALL_obs_p4.sess = cat(1,ALL_obs_p4.sess,iFD*ones(size(temp4)));
   
   % print something about this session   
   %fprintf('%s (%s): p0 L (N = %d) mean %.2f +/- %.2f, median %.2f; R (N = %d) mean %.2f +/- %.2f, median %.2f\n', ...
       %ExpKeys.goodSWR{1}(1:15),cfg.sessions{this_session_type},nL,nanmean(tempL0),nanstd(tempL0),nanmedian(tempL0),nR,nanmean(tempR0),nanstd(tempR0),nanmedian(tempR0));
   
   %fprintf('%s (%s): p4 L (N = %d) mean %.2f +/- %.2f, median %.2f; R (N = %d) mean %.2f +/- %.2f, median %.2f\n', ...
       %ExpKeys.goodSWR{1}(1:15),cfg.sessions{this_session_type},nL,nanmean(tempL4),nanstd(tempL4),nanmedian(tempL4),nR,nanmean(tempR4),nanstd(tempR4),nanmedian(tempR4));
   
   
end

%% plot observed data (p0)
p0_food_left = ALL_obs_p0.data(ALL_obs_p0.arm == 1 & ALL_obs_p0.type == 1);
p0_food_right = ALL_obs_p0.data(ALL_obs_p0.arm == 2 & ALL_obs_p0.type == 1);

p0_water_left = ALL_obs_p0.data(ALL_obs_p0.arm == 1 & ALL_obs_p0.type == 2);
p0_water_right = ALL_obs_p0.data(ALL_obs_p0.arm == 2 & ALL_obs_p0.type == 2);

data.all.p0_food_left = p0_food_left;
data.all.p0_food_right = p0_food_right;
data.all.p0_water_left = p0_water_left;
data.all.p0_water_right = p0_water_right;

disp(' ') 
disp('Displaying mean p0:')
disp(['Left singles, food restricted:', num2str(nanmean(p0_food_left))])
disp(['Right singles, food restricted:', num2str(nanmean(p0_food_right))])
disp(['Left singles, water restricted:', num2str(nanmean(p0_water_left))])
disp(['Right singles, water restricted:', num2str(nanmean(p0_water_right))])
disp(' ')

% figure; hold on
% 
% xbin = [1 2 4 5];
% d = [nanmean(p0_food_left) nanmean(p0_food_right) nanmean(p0_water_left) nanmean(p0_water_right)];
% col = {chooseFood_color chooseFood_color chooseWater_color chooseWater_color };
% 
% for iB = 1:length(xbin)
%    
%     %h(iB) = bar(xbin(iB),d(iB)); set(h(iB),'FaceColor',col{iB},'EdgeColor','none');
%     hold on;
%     
% end

% if min(d) < 0
%     ymin = min(d) - 0.03;
% else
%     ymin = 0;
% end
% 
% if max(d) > 0.3
%     ymax = max(d) + 0.03;
% else
%     ymax = 0.3;
% end
% 
% set(gca,'XTick',xbin,'XTickLabel',{'left','right','left','right'}, ...
%     'YLim',[ymin ymax],'LineWidth',1,'YTick',0:0.05:0.3,'FontSize',14);
% 
% 
% 
% ylabel('proportion of HFEs active');
% xlabel('  (food restricted)                                 (water restricted)')
% title(sprintf('%d rats, %d sessions, %d left cells, %d right cells',length(cfg.rats),length(fd),length(p0_food_left)+length(p0_water_left),length(p0_water_right)+length(p0_food_right)));
% 
% cd(cfg.output_fd)
% if cfg.printFig
%     print(gcf,'-dpng','-r300',cat(2,'Collect_CoOccur_p0',cfg.output_suffix,'.png'));
% end

%% p0 stats
%statistics regarding the fraction of Q-matrix time bins each cell was active in
% more specifically, the fraction of candidate events each cell was active
% in alone

% food restrict: any difference left vs right?
d1 = ALL_obs_p0.data(ALL_obs_p0.type == 1 & ALL_obs_p0.arm == 1);
d2 = ALL_obs_p0.data(ALL_obs_p0.type == 1 & ALL_obs_p0.arm == 2);
fprintf('left v right (food): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
    nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));

% water restrict: any difference left vs right?
d1 = ALL_obs_p0.data(ALL_obs_p0.type == 2 & ALL_obs_p0.arm == 1);
d2 = ALL_obs_p0.data(ALL_obs_p0.type == 2 & ALL_obs_p0.arm == 2);
fprintf('left v right (water): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
    nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));

% interaction
Y = ALL_obs_p0.data; % this is p-value data about the fraction of SWRs each cell was active in by itself
X1 = ALL_obs_p0.arm; % this identifies whether the p-val belongs to a left cell or a right cell
X2 = ALL_obs_p0.type; % this identifies whether the cell was recorded on a food-restriction or a water-restriction day

[P,T,STATS,TERMS] = anovan(Y,cat(2,X1,X2),'model','full','varnames',{'arm','restr'},'continuous',[]);

data.all.anova.P = P;

disp(' ')
disp('*******************************************************************')
disp('***                                                             ***')
disp('***   Displaying ANOVA results for single cell activity (p0):   ***')
disp('***                                                             ***')
disp('*******************************************************************')
disp(' ')
disp(T)

%% plot observed data (p4)
p4_food_left = ALL_obs_p4.data(ALL_obs_p4.arm == 1 & ALL_obs_p4.type == 1);
p4_food_right = ALL_obs_p4.data(ALL_obs_p4.arm == 2 & ALL_obs_p4.type == 1);

p4_water_left = ALL_obs_p4.data(ALL_obs_p4.arm == 1 & ALL_obs_p4.type == 2);
p4_water_right = ALL_obs_p4.data(ALL_obs_p4.arm == 2 & ALL_obs_p4.type == 2);

data.all.p4_food_left = p4_food_left;
data.all.p4_food_right = p4_food_right;
data.all.p4_water_left = p4_water_left;
data.all.p4_water_right = p4_water_right;

disp(' ') 
disp('Displaying mean p4:')
disp(['Left pairs, food restricted:', num2str(nanmean(p4_food_left))])
disp(['Right pairs, food restricted:', num2str(nanmean(p4_food_right))])
disp(['Left pairs, water restricted:', num2str(nanmean(p4_water_left))])
disp(['Right pairs, water restricted:', num2str(nanmean(p4_water_right))])
disp(' ')

% figure; hold on
% 
% xbin = [1 2 4 5];
d = [nanmean(p4_food_left) nanmean(p4_food_right) nanmean(p4_water_left) nanmean(p4_water_right)];
% col = {chooseFood_color chooseFood_color chooseWater_color chooseWater_color };
% 
% for iB = 1:length(xbin)
%    
%     %h(iB) = bar(xbin(iB),d(iB)); set(h(iB),'FaceColor',col{iB},'EdgeColor','none');
%     hold on;
%     
% end
% 
% if min(d) < 0
%     ymin = min(d) - 0.03;
% else
%     ymin = 0;
% end
% 
% if max(d) >= 0.95
%     ymax = 1.03;
% else
%     ymax = 1;
% end
% 
% set(gca,'XTick',xbin,'XTickLabel',{'left','right','left','right'}, ...
%     'YLim',[ymin ymax],'LineWidth',1,'YTick',0:0.25:1,'FontSize',14);
% 
% ylabel('HFE coactivation Z-score');
% xlabel('  (food restricted)                                 (water restricted)')
% title(sprintf('%d rats, %d sessions, %d left pairs, %d right pairs',length(cfg.rats), length(fd), length(p4_food_left)+length(p4_water_left),length(p4_food_right)+length(p4_food_right)));
% 
% if cfg.printFig
%     print(gcf,'-dpng','-r300',cat(2,'Collect_CoOccur_p4',cfg.output_suffix,'.png'));
% end

%% p4 stats

% food restrict: any difference left vs right?
d1 = ALL_obs_p4.data(ALL_obs_p4.type == 1 & ALL_obs_p4.arm == 1);
d2 = ALL_obs_p4.data(ALL_obs_p4.type == 1 & ALL_obs_p4.arm == 2);
fprintf('left v right (food): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
    nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));

% water restrict: any difference left vs right?
d1 = ALL_obs_p4.data(ALL_obs_p4.type == 2 & ALL_obs_p4.arm == 1);
d2 = ALL_obs_p4.data(ALL_obs_p4.type == 2 & ALL_obs_p4.arm == 2);
fprintf('left v right (water): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
    nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));

% interaction
Y = ALL_obs_p4.data; % z-score of observed coactivity of each cell pair as compared to randomly shuffled data 
X1 = ALL_obs_p4.arm; % this identifies whether the p-val belongs to a left pair or a right pair
X2 = ALL_obs_p4.type; % this identifies whether the pair was recorded on a food-restriction or a water-restriction day


[P,T,STATS,TERMS] = anovan(Y,cat(2,X1,X2),'model','full','varnames',{'arm','restr'},'continuous',[]);
disp(' ')
disp('*******************************************************************')
disp('***                                                             ***')
disp('***   Displaying ANOVA results for cell pair coactivity (p4):   ***')
disp('***                                                             ***')
disp('*******************************************************************')
disp(' ')
disp(T)

%% GET INDIVIDUAL DATA

for iRat = 1:length(cfg.rats)
    disp(' ')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%                                                               %%%')
    disp(['%%%                Co-occurrence data for ',cfg.rats{iRat},':                   %%%'])
    disp('%%%                                                               %%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); disp(' ')
    cfg_temp.rats = cfg.rats(iRat);
    fd = getTmazeDataPath(cfg_temp);
    
    ALL_obs_p0.data = []; % the actual observed activation probabilities
    ALL_obs_p0.arm = []; % left (1), right (2)
    ALL_obs_p0.type = []; % restriction type (food/water)
    ALL_obs_p0.sess = []; % session ID in case we want to restrict later
    
    ALL_obs_p4.data = []; % the actual observed activation probabilities
    ALL_obs_p4.arm = []; % left (1), right (2)
    ALL_obs_p4.type = []; % restriction type (food/water)
    ALL_obs_p4.sess = []; % session ID in case we want to restrict later
    
    for iFD = 1:length(fd)
        
        cd(fd{iFD});
        LoadExpKeys;
        
        cd([pwd,'\','files']);
        this_file = FindFile(cat(2,cfg.prefix,'*coOccur_data.mat'));
        load(this_file);
        
        %
        this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
        
        % collect data
        
        % get data to process from loaded file
        this_data = out.(cfg.what); %eval(cat(2,'out.',cfg.what,';'));
        nL = size(this_data.left,1); % number of left cells
        nR = size(this_data.right,1); % number of right cells
        
        tempL0 = this_data.full_cc.p0(1:nL); tempR0 = this_data.full_cc.p0(nL+1:end);
        if size(tempL0,2) == 0, tempL0 = tempL0'; end; if size(tempR0,2) == 0, tempR0 = tempR0'; end
        temp0 = cat(1,tempL0,tempR0);
        ALL_obs_p0.data = cat(1,ALL_obs_p0.data,temp0); % the actual observed activation probabilities
        ALL_obs_p0.arm = cat(1,ALL_obs_p0.arm,cat(1,ones(length(tempL0),1),2*ones(length(tempR0),1))); % this make an identity array
        % thing that contains 1's in places corresponding to left cell data (probabilities), and 2's for right cell data
        ALL_obs_p0.type = cat(1,ALL_obs_p0.type,this_session_type*ones(size(temp0)));
        ALL_obs_p0.sess = cat(1,ALL_obs_p0.sess,iFD*ones(size(temp0)));
        
        tempL4 = this_data.left(:); tempR4 = this_data.right(:); temp4 = cat(1,tempL4,tempR4);
        ALL_obs_p4.data = cat(1,ALL_obs_p4.data,temp4); % the actual observed activation probabilities
        ALL_obs_p4.arm = cat(1,ALL_obs_p4.arm,cat(1,ones(length(tempL4),1),2*ones(length(tempR4),1)));
        ALL_obs_p4.type = cat(1,ALL_obs_p4.type,this_session_type*ones(size(temp4)));
        ALL_obs_p4.sess = cat(1,ALL_obs_p4.sess,iFD*ones(size(temp4)));
        
        % print something about this session
        fprintf('%s (%s): p0 L (N = %d) mean %.2f +/- %.2f, median %.2f; R (N = %d) mean %.2f +/- %.2f, median %.2f\n', ...
            ExpKeys.goodSWR{1}(1:15),cfg.sessions{this_session_type},nL,nanmean(tempL0),nanstd(tempL0),nanmedian(tempL0),nR,nanmean(tempR0),nanstd(tempR0),nanmedian(tempR0));
        
        fprintf('%s (%s): p4 L (N = %d) mean %.2f +/- %.2f, median %.2f; R (N = %d) mean %.2f +/- %.2f, median %.2f\n', ...
            ExpKeys.goodSWR{1}(1:15),cfg.sessions{this_session_type},nL,nanmean(tempL4),nanstd(tempL4),nanmedian(tempL4),nR,nanmean(tempR4),nanstd(tempR4),nanmedian(tempR4));
        
        
    end
    
    %% plot observed data (p0)
    p0_food_left = ALL_obs_p0.data(ALL_obs_p0.arm == 1 & ALL_obs_p0.type == 1);
    p0_food_right = ALL_obs_p0.data(ALL_obs_p0.arm == 2 & ALL_obs_p0.type == 1);
    
    p0_water_left = ALL_obs_p0.data(ALL_obs_p0.arm == 1 & ALL_obs_p0.type == 2);
    p0_water_right = ALL_obs_p0.data(ALL_obs_p0.arm == 2 & ALL_obs_p0.type == 2);
    
    data.(cfg.rats{iRat}).p0_food_left = p0_food_left;
    data.(cfg.rats{iRat}).p0_food_right = p0_food_right;
    data.(cfg.rats{iRat}).p0_water_left = p0_water_left;
    data.(cfg.rats{iRat}).p0_water_right = p0_water_right;
    
    disp(' ')
    disp('Displaying mean p0:')
    disp(['Left singles, food restricted:', num2str(nanmean(p0_food_left))])
    disp(['Right singles, food restricted:', num2str(nanmean(p0_food_right))])
    disp(['Left singles, water restricted:', num2str(nanmean(p0_water_left))])
    disp(['Right singles, water restricted:', num2str(nanmean(p0_water_right))])
    disp(' ')
    
%     figure; hold on
%     
%     xbin = [1 2 4 5];
    d = [nanmean(p0_food_left) nanmean(p0_food_right) nanmean(p0_water_left) nanmean(p0_water_right)];
%     col = {chooseFood_color chooseFood_color chooseWater_color chooseWater_color };
%     
%     for iB = 1:length(xbin)
%         
%         %h(iB) = bar(xbin(iB),d(iB)); set(h(iB),'FaceColor',col{iB},'EdgeColor','none');
%         hold on;
%         
%     end
%     
%     if min(d) < 0
%         ymin = min(d) - 0.03;
%     else
%         ymin = 0;
%     end
%     
%     if max(d) > 0.3
%         ymax = max(d) + 0.03;
%     else
%         ymax = 0.3;
%     end
%     
%     set(gca,'XTick',xbin,'XTickLabel',{'left','right','left','right'}, ...
%         'YLim',[ymin ymax],'LineWidth',1,'YTick',0:0.05:0.3,'FontSize',14);
%     
%     
%     
%     ylabel('proportion of HFEs active');
%     xlabel('  (food restricted)                                 (water restricted)')
%     title(sprintf('%d rats, %d sessions, %d left cells, %d right cells',length(cfg.rats),length(fd),length(p0_food_left)+length(p0_water_left),length(p0_water_right)+length(p0_food_right)));
%     
%     cd(cfg.output_fd)
%     if cfg.printFig
%         print(gcf,'-dpng','-r300',cat(2,'Collect_CoOccur_p0',cfg.output_suffix,'.png'));
%     end
    
    %% p0 stats
    %statistics regarding the fraction of Q-matrix time bins each cell was active in
    % more specifically, the fraction of candidate events each cell was active
    % in alone
    
    % food restrict: any difference left vs right?
    d1 = ALL_obs_p0.data(ALL_obs_p0.type == 1 & ALL_obs_p0.arm == 1);
    d2 = ALL_obs_p0.data(ALL_obs_p0.type == 1 & ALL_obs_p0.arm == 2);
    fprintf('left v right (food): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
        nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));
    
    % water restrict: any difference left vs right?
    d1 = ALL_obs_p0.data(ALL_obs_p0.type == 2 & ALL_obs_p0.arm == 1);
    d2 = ALL_obs_p0.data(ALL_obs_p0.type == 2 & ALL_obs_p0.arm == 2);
    fprintf('left v right (water): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
        nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));
    
    % interaction
    Y = ALL_obs_p0.data; % this is p-value data about the fraction of SWRs each cell was active in by itself
    X1 = ALL_obs_p0.arm; % this identifies whether the p-val belongs to a left cell or a right cell
    X2 = ALL_obs_p0.type; % this identifies whether the cell was recorded on a food-restriction or a water-restriction day
    
    [P,T,STATS,TERMS] = anovan(Y,cat(2,X1,X2),'model','full','varnames',{'arm','restr'},'continuous',[]);
    
    data.all.anova.P = P;
    
    disp(' ')
    disp('*******************************************************************')
    disp('***                                                             ***')
    disp('***   Displaying ANOVA results for single cell activity (p0):   ***')
    disp('***                                                             ***')
    disp('*******************************************************************')
    disp(' ')
    disp(T)
    
    %% plot observed data (p4)
    p4_food_left = ALL_obs_p4.data(ALL_obs_p4.arm == 1 & ALL_obs_p4.type == 1);
    p4_food_right = ALL_obs_p4.data(ALL_obs_p4.arm == 2 & ALL_obs_p4.type == 1);
    
    p4_water_left = ALL_obs_p4.data(ALL_obs_p4.arm == 1 & ALL_obs_p4.type == 2);
    p4_water_right = ALL_obs_p4.data(ALL_obs_p4.arm == 2 & ALL_obs_p4.type == 2);
    
    data.(cfg.rats{iRat}).p4_food_left = p4_food_left;
    data.(cfg.rats{iRat}).p4_food_right = p4_food_right;
    data.(cfg.rats{iRat}).p4_water_left = p4_water_left;
    data.(cfg.rats{iRat}).p4_water_right = p4_water_right;
    
    disp(' ')
    disp('Displaying mean p4:')
    disp(['Left pairs, food restricted:', num2str(nanmean(p4_food_left))])
    disp(['Right pairs, food restricted:', num2str(nanmean(p4_food_right))])
    disp(['Left pairs, water restricted:', num2str(nanmean(p4_water_left))])
    disp(['Right pairs, water restricted:', num2str(nanmean(p4_water_right))])
    disp(' ')
    
%     figure; hold on
%     
%     xbin = [1 2 4 5];
    d = [nanmean(p4_food_left) nanmean(p4_food_right) nanmean(p4_water_left) nanmean(p4_water_right)];
%     col = {chooseFood_color chooseFood_color chooseWater_color chooseWater_color };
%     
%     for iB = 1:length(xbin)
%         
%         %h(iB) = bar(xbin(iB),d(iB)); set(h(iB),'FaceColor',col{iB},'EdgeColor','none');
%         hold on;
%         
%     end
%     
%     if min(d) < 0
%         ymin = min(d) - 0.03;
%     else
%         ymin = 0;
%     end
%     
%     if max(d) >= 0.95
%         ymax = 1.03;
%     else
%         ymax = 1;
%     end
%     
%     set(gca,'XTick',xbin,'XTickLabel',{'left','right','left','right'}, ...
%         'YLim',[ymin ymax],'LineWidth',1,'YTick',0:0.25:1,'FontSize',14);
%     
%     ylabel('HFE coactivation Z-score');
%     xlabel('  (food restricted)                                 (water restricted)')
%     title(sprintf('%d rats, %d sessions, %d left pairs, %d right pairs',length(cfg.rats), length(fd), length(p4_food_left)+length(p4_water_left),length(p4_food_right)+length(p4_food_right)));
%     
%     if cfg.printFig
%         print(gcf,'-dpng','-r300',cat(2,'Collect_CoOccur_p4',cfg.output_suffix,'.png'));
%     end
    
    %% p4 stats
    
    % food restrict: any difference left vs right?
    d1 = ALL_obs_p4.data(ALL_obs_p4.type == 1 & ALL_obs_p4.arm == 1);
    d2 = ALL_obs_p4.data(ALL_obs_p4.type == 1 & ALL_obs_p4.arm == 2);
    fprintf('left v right (food): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
        nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));
    
    % water restrict: any difference left vs right?
    d1 = ALL_obs_p4.data(ALL_obs_p4.type == 2 & ALL_obs_p4.arm == 1);
    d2 = ALL_obs_p4.data(ALL_obs_p4.type == 2 & ALL_obs_p4.arm == 2);
    fprintf('left v right (water): %.2f +/- %.2f v %.2f +/- %.2f, p(%d,%d) = %.2e\n', ...
        nanmean(d1),nanstd(d1),nanmean(d2),nanstd(d2),length(d1),length(d2),ranksum(d1,d2));
    
    % interaction
    Y = ALL_obs_p4.data; % z-score of observed coactivity of each cell pair as compared to randomly shuffled data
    X1 = ALL_obs_p4.arm; % this identifies whether the p-val belongs to a left pair or a right pair
    X2 = ALL_obs_p4.type; % this identifies whether the pair was recorded on a food-restriction or a water-restriction day
    
    
    [P,T,STATS,TERMS] = anovan(Y,cat(2,X1,X2),'model','full','varnames',{'arm','restr'},'continuous',[]);
    disp(' ')
    disp('*******************************************************************')
    disp('***                                                             ***')
    disp('***   Displaying ANOVA results for cell pair coactivity (p4):   ***')
    disp('***                                                             ***')
    disp('*******************************************************************')
    disp(' ')
    disp(T)
    
end
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

if cfg.writeDiary, diary off, end

if cfg.saveData; cd(cfg.output_fd); save([cfg.prefix,'cooccurrence_rush.mat'],'data'); end
cd(originalFolder)
