% batch script to collect and plot co-occurrence data
% SfN poster rush code
% aacarey Oct 2015, modified from MvdM's ALL_CollectCorrSequences

clear

%%

warning off

cfg.prefix = 'CS_ALL_matched_'; % which files to load

% which events to process?
cfg.whichEvents = 'all';% 'all','prerecord','task','postrecord'

cfg.p = 0.99; % p-level for which events to count
cfg.sessions = {'food','water'};
cfg.arms = {'left','right'}; % needs to match order of out.score1[] and out.score3[] (1 is left, 2 right arm)

cfg.writeDiary = 1; % keep a text file record of command window history

cfg.output_fd = 'E:\Documents\TmazePaper\data';

% do you want to save the data?
cfg.saveData = 1;

%%
fg = []; cfg.rats = {'R042','R044','R050','R064'};
cfg.requireCandidates = 1;
fd = getTmazeDataPath(cfg);

iWasHere = pwd;
switch cfg.whichEvents
    case 'all'
        cfg.output_prefix = 'CS_ALL_matched_';
    case 'prerecord'
        cfg.output_prefix = 'CS_PRE_matched_';
    case 'postrecord'
        cfg.output_prefix = 'CS_POST_matched_';
        
    case 'taskrest'
        cfg.output_prefix = 'CS_TASKREST_';
        
    case 'task'
        cfg.output_prefix = 'CS_TASK_matched_';
    otherwise
        error('Unrecognized cfg.whichEvents')
end


if cfg.writeDiary % save command window text
    warning off
    cd(cfg.output_fd)
    diary([cfg.output_prefix,'stats_rush','.txt'])
    cd(cfg.output_fd)
    disp(' ')
    disp(date)
    disp(' ')
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
disp('%%%               Combined sequences data (all rats):             %%%')
disp('%%%                                                               %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% init vars to place data into

ALL_sig_seq.count = []; % the actual observed activation probabilities
ALL_sig_seq.arm = []; % left (1), right (2)
ALL_sig_seq.type = []; % restriction type (food/water)
ALL_sig_seq.sess = []; % session ID in case we want to restrict later

%% collect
nSessions = 0;
for iFD = 1:length(fd)
    nSessions = nSessions +1;
    close all;
    
   cd(fd{iFD}); 
   LoadExpKeys;
   
   %cd('files');
   cd([fd{iFD},'\files'])
   [~,sessionID,~] = fileparts(fd{iFD}); % 'RXXX-201X-XX-XX'
   this_file = FindFiles([cfg.prefix,sessionID,'-CorrScores.mat']);
   
   if isempty(this_file)
      fprintf('Session %s: no CorrScores file found, skipping...\n',fd{iFD}); 
      continue;
   end
   
   load(this_file{1}); cd ..
   
   this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
   
   % simply count the number of significant sequences
   for iArm = 1:length(cfg.arms)
       
       % optionally, restrict to specific task epoch
       
       keep_idx{iArm} = find(out.score1(iArm).WIN_rho_perc > cfg.p & out.score3(iArm).WIN_rho_perc > cfg.p);
   
       switch cfg.whichEvents
           case 'prerecord'
               t = [out.score1(iArm).WIN_iv(:).tstart];
               t = t(keep_idx{iArm});
               keep_idx{iArm} = keep_idx{iArm}(t < ExpKeys.prerecord(2));
           case 'task'
               t = [out.score1(iArm).WIN_iv(:).tstart];
               t = t(keep_idx{iArm});
               keep_idx{iArm} = keep_idx{iArm}(t > ExpKeys.task(1) & t < ExpKeys.task(2));
           case 'postrecord'
               t = [out.score1(iArm).WIN_iv(:).tstart];
               t = t(keep_idx{iArm});
               keep_idx{iArm} = keep_idx{iArm}(t > ExpKeys.postrecord(1));
       end
                 
       ALL_sig_seq.count = cat(1,ALL_sig_seq.count,length(keep_idx{iArm}));
       ALL_sig_seq.arm = cat(1,ALL_sig_seq.arm,iArm);
       ALL_sig_seq.type = cat(1,ALL_sig_seq.type,this_session_type);
       ALL_sig_seq.sess = cat(1,ALL_sig_seq.sess,iFD);
       
   end
   
   % output something about this session
   %fprintf('%s (%s): p0 L %d, R %d\n', ...
      % ExpKeys.goodSWR{1}(1:15),cfg.sessions{this_session_type},length(keep_idx{1}),length(keep_idx{2}));  
end
 

food_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 1);
food_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 1);

water_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 2);
water_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 2);

d = [nansum(food_left) nansum(food_right) nansum(water_left) nansum(water_right)];

data.all.food_left = nansum(food_left);
data.all.food_right = nansum(food_right);
data.all.water_left = nansum(water_left);
data.all.water_right = nansum(water_right);
disp(' ')
disp('Displaying number of significant sequences: ')
disp(['Left, food restr, ',num2str(data.all.food_left)])
disp(['Right, food restr ',num2str(data.all.food_right)])   
disp(['Left, water restr ',num2str(data.all.water_left)])
disp(['Right, water restr ',num2str(data.all.water_right)])
disp(' ')

p_foodseq = tmaze_coin_pval(max([d(1) d(2)]),d(1)+d(2));
p_waterseq = tmaze_coin_pval(max([d(3) d(4)]),d(3)+d(4));
disp(' '); 
disp('***************************')
disp('** binomial test results **')
disp('***************************')
disp(['Food day pval: ',num2str(p_foodseq)])
disp(['Water day pval: ',num2str(p_waterseq)])

temp = data.all; % because I don't want to type as much down there
    
    % [nleft on food day, nright on food day, nleft on water day, nright on water day]
    obs = [temp.food_left temp.food_right temp.water_left temp.water_right]; % the observed, real data
    
    nFood = obs(1) + obs(2); % total number of food day trials. could also do temp.nFood
    nWater = obs(3) + obs(4); % total number of water day trials
    ratioL = (obs(1)+obs(3)) / sum(obs); % ratio of L trials compared to total trials
    ratioR = (obs(2)+obs(4)) / sum(obs); % ratio of R trials
    
    exp = [nFood*ratioL nFood*ratioR nWater*ratioL nWater*ratioR]; % the "expected" values
    
    bins = 0:3;
    edges = -0.5:3.5;
    disp(' '); 
    disp('***************************************')
    disp('**                                   **')
    disp('**      chi square test results      **')
    disp('**                                   **')
    disp('***************************************')
    [h,p,stats] = chi2gof(bins,'Alpha',0.01,'Edges',edges,'freq',obs,'expected',exp,'Emin',1) % no semicolon b/c want to see output

%% GET INDIVIDUAL DATA
rats = cfg.rats; 
for iRat = 1:length(cfg.rats)
    disp(' ')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%                                                               %%%')
    disp(['%%%                  Sequences data for ',cfg.rats{iRat},':                     %%%'])
    disp('%%%                                                               %%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'); disp(' ')
    cfg_temp.rats = cfg.rats(iRat);
    fd = getTmazeDataPath(cfg_temp);
    
    %% init vars to place data into
    
    ALL_sig_seq.count = []; % the actual observed activation probabilities
    ALL_sig_seq.arm = []; % left (1), right (2)
    ALL_sig_seq.type = []; % restriction type (food/water)
    ALL_sig_seq.sess = []; % session ID in case we want to restrict later
    
    %% collect
    nSessions = 0;
    for iFD = 1:length(fd)
        disp(' ')
        nSessions = nSessions +1;
        close all;
        
        cd(fd{iFD});
        LoadExpKeys;
        
        %cd('files');
        cd([fd{iFD},'\files'])
        [~,sessionID,~] = fileparts(fd{iFD}); % 'RXXX-201X-XX-XX'
        this_file = FindFiles([cfg.prefix,sessionID,'-CorrScores.mat']);
        
        if isempty(this_file)
            fprintf('Session %s: no CorrScores file found, skipping...\n',fd{iFD});
            continue;
        end
        
        load(this_file{1}); cd ..
        
        this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
        
        % simply count the number of significant sequences
        for iArm = 1:length(cfg.arms)
            
            % optionally, restrict to specific task epoch
            
            keep_idx{iArm} = find(out.score1(iArm).WIN_rho_perc > cfg.p & out.score3(iArm).WIN_rho_perc > cfg.p);
            
            switch cfg.whichEvents
                case 'prerecord'
                    t = [out.score1(iArm).WIN_iv(:).tstart];
                    t = t(keep_idx{iArm});
                    keep_idx{iArm} = keep_idx{iArm}(t < ExpKeys.prerecord(2));
                case 'task'
                    t = [out.score1(iArm).WIN_iv(:).tstart];
                    t = t(keep_idx{iArm});
                    keep_idx{iArm} = keep_idx{iArm}(t > ExpKeys.task(1) & t < ExpKeys.task(2));
                case 'postrecord'
                    t = [out.score1(iArm).WIN_iv(:).tstart];
                    t = t(keep_idx{iArm});
                    keep_idx{iArm} = keep_idx{iArm}(t > ExpKeys.postrecord(1));
            end
            
            ALL_sig_seq.count = cat(1,ALL_sig_seq.count,length(keep_idx{iArm}));
            ALL_sig_seq.arm = cat(1,ALL_sig_seq.arm,iArm);
            ALL_sig_seq.type = cat(1,ALL_sig_seq.type,this_session_type);
            ALL_sig_seq.sess = cat(1,ALL_sig_seq.sess,iFD);
            
        end
        
        % output something about this session
        fprintf('%s (%s): p0 L %d, R %d\n', ...
            ExpKeys.goodSWR{1}(1:15),cfg.sessions{this_session_type},length(keep_idx{1}),length(keep_idx{2}));
    end
    
    
    food_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 1);
    food_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 1);
    
    water_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 2);
    water_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 2);
    
    d = [nansum(food_left) nansum(food_right) nansum(water_left) nansum(water_right)];
    
    data.(rats{iRat}).food_left = nansum(food_left);
    data.(rats{iRat}).food_right = nansum(food_right);
    data.(rats{iRat}).water_left = nansum(water_left);
    data.(rats{iRat}).water_right = nansum(water_right);
    
    disp(' ')
    disp('Displaying number of significant sequences: ')
    disp(['Left, food restr, ',num2str(data.(rats{iRat}).food_left)])
    disp(['Right, food restr ',num2str(data.(rats{iRat}).food_right)])
    disp(['Left, water restr ',num2str(data.(rats{iRat}).water_left)])
    disp(['Right, water restr ',num2str(data.(rats{iRat}).water_right)])
    disp(' ')
    
    p_foodseq = tmaze_coin_pval(max([d(1) d(2)]),d(1)+d(2));
    p_waterseq = tmaze_coin_pval(max([d(3) d(4)]),d(3)+d(4));
    disp(' ');
    disp('***************************')
    disp('** binomial test results **')
    disp('***************************')
    disp(['Food day pval: ',num2str(p_foodseq)])
    disp(['Water day pval: ',num2str(p_waterseq)])
    
   
    temp = data.(rats{iRat}); % because I don't want to type as much down there
    
    % [nleft on food day, nright on food day, nleft on water day, nright on water day]
    obs = [temp.food_left temp.food_right temp.water_left temp.water_right]; % the observed, real data
    
    nFood = obs(1) + obs(2); % total number of food day trials. could also do temp.nFood
    nWater = obs(3) + obs(4); % total number of water day trials
    ratioL = (obs(1)+obs(3)) / sum(obs); % ratio of L trials compared to total trials
    ratioR = (obs(2)+obs(4)) / sum(obs); % ratio of R trials
    
    exp = [nFood*ratioL nFood*ratioR nWater*ratioL nWater*ratioR]; % the "expected" values
    
    bins = 0:3;
    edges = -0.5:3.5;
    disp(' ');
    disp('***************************************')
    disp('**                                   **')
    disp('**      chi square test results      **')
    disp('**                                   **')
    disp('***************************************')
    [h,p,stats] = chi2gof(bins,'Alpha',0.01,'Edges',edges,'freq',obs,'expected',exp,'Emin',1) % no semicolon b/c want to see output
end

data.date = datestr(now);
data.script = mfilename;

if cfg.writeDiary; diary off; end

if cfg.saveData; cd(cfg.output_fd); save([cfg.output_prefix,'sequences_rush.mat'],'data'); end

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

cd(iWasHere)
warning on