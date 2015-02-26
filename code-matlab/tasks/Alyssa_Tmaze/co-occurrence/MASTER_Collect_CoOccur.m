% batch script to collect and plot co-occurrence data

%%
cfg = []; cfg.rats = {'R042','R044','R050'};
cfg.requireCandidates = 1;
fd = getTmazeDataPath(cfg);

cfg.prefix = 'CC_PRE';
cfg.what = 'SWRz';

cfg.sessions = {'food','water'};
cfg.arms = {'left','right'}; % to make ALL_obs_p0.arm

%% init vars to place data into

ALL_obs_p0.data = []; % the actual observed activation probabiliries
ALL_obs_p0.arm = []; % left (1), right (2)
ALL_obs_p0.type = []; % restriction type (food/water)
ALL_obs_p0.sess = []; % session ID in case we want to restrict later

ALL_obs_p4.data = []; % the actual observed activation probabiliries
ALL_obs_p4.arm = []; % left (1), right (2)
ALL_obs_p4.type = []; % restriction type (food/water)
ALL_obs_p4.sess = []; % session ID in case we want to restrict later


%%
for iFD = 1:length(fd)
    
   cd(fd{iFD}); 
   LoadExpKeys;
   
   cd('files');
   this_file = FindFile(cat(2,cfg.prefix,'*coOccur_data.mat'));
   load(this_file);
   
   %
   this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));

   % collect data
   
   % get data to process from loaded file
   this_data = eval(cat(2,'out.',cfg.what,';'));
   nL = size(this_data.left,1); nR = size(this_data.right,1);
   
   tempL0 = this_data.full_cc.p0(1:nL); tempR0 = this_data.full_cc.p0(nL+1:end); 
   if size(tempL0,2) == 0, tempL0 = tempL0'; end; if size(tempR0,2) == 0, tempR0 = tempR0'; end
   temp0 = cat(1,tempL0,tempR0);
   ALL_obs_p0.data = cat(1,ALL_obs_p0.data,temp0); % the actual observed activation probabiliries
   ALL_obs_p0.arm = cat(1,ALL_obs_p0.arm,cat(1,ones(length(tempL0),1),2*ones(length(tempR0),1)));
   ALL_obs_p0.type = cat(1,ALL_obs_p0.type,this_session_type*ones(size(temp0)));
   ALL_obs_p0.sess = cat(1,ALL_obs_p0.sess,iFD*ones(size(temp0)));
   
   tempL4 = this_data.left(:); tempR4 = this_data.right(:); temp4 = cat(1,tempL4,tempR4);
   ALL_obs_p4.data = cat(1,ALL_obs_p4.data,temp4); % the actual observed activation probabiliries
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

subplot(221);

xbin = [1 2 4 5];
d = [nanmean(p0_food_left) nanmean(p0_food_right) nanmean(p0_water_left) nanmean(p0_water_right)];
cl = 'rrbb';

for iB = 1:length(xbin)
   
    h(iB) = bar(xbin(iB),d(iB)); set(h(iB),'FaceColor',cl(iB),'EdgeColor','none');
    hold on;
    
end

set(gca,'XTick',xbin,'XTickLabel',{'left','right','left','right'}, ...
    'YLim',[0 0.3],'LineWidth',1,'YTick',0:0.05:0.3,'FontSize',18);

ylabel('proportion of SWRs active');
xlabel('(food restricted)         (water restricted)');
title(sprintf('%d rats, %d sessions, %d left cells, %d right cells',length(cfg.rats),length(fd),length(p0_food_left)+length(p0_water_left),length(p0_food_right)+length(p0_food_right)));

%% p0 stats

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
Y = ALL_obs_p0.data;
X1 = ALL_obs_p0.arm; % left/right
X2 = ALL_obs_p0.type; % food/water

[P,T,STATS,TERMS] = anovan(Y,cat(2,X1,X2),'model','full','varnames',{'arm','restr'},'continuous',[]);

%% plot observed data (p4)
p4_food_left = ALL_obs_p4.data(ALL_obs_p4.arm == 1 & ALL_obs_p4.type == 1);
p4_food_right = ALL_obs_p4.data(ALL_obs_p4.arm == 2 & ALL_obs_p4.type == 1);

p4_water_left = ALL_obs_p4.data(ALL_obs_p4.arm == 1 & ALL_obs_p4.type == 2);
p4_water_right = ALL_obs_p4.data(ALL_obs_p4.arm == 2 & ALL_obs_p4.type == 2);

subplot(222);

xbin = [1 2 4 5];
d = [nanmean(p4_food_left) nanmean(p4_food_right) nanmean(p4_water_left) nanmean(p4_water_right)];
cl = 'rrbb';

for iB = 1:length(xbin)
   
    h(iB) = bar(xbin(iB),d(iB)); set(h(iB),'FaceColor',cl(iB),'EdgeColor','none');
    hold on;
    
end

set(gca,'XTick',xbin,'XTickLabel',{'left','right','left','right'}, ...
    'YLim',[0 1],'LineWidth',1,'YTick',0:0.25:1,'FontSize',18);

ylabel('SWR coactivation Z-score');
xlabel('(food restricted)         (water restricted)');
title(sprintf('%d rats, %d sessions, %d left pairs, %d right pairs',length(cfg.rats),length(fd),length(p4_food_left)+length(p4_water_left),length(p4_food_right)+length(p4_food_right)));

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
Y = ALL_obs_p4.data;
X1 = ALL_obs_p4.arm; % left/right
X2 = ALL_obs_p4.type; % food/water

[P,T,STATS,TERMS] = anovan(Y,cat(2,X1,X2),'model','full','varnames',{'arm','restr'},'continuous',[]);
