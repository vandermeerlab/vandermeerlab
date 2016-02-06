%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                   Collect Arm Choice Behaviour                      %%%                
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Tmaze data

% ABOUT
% This script loads metadata and ExpKeys and collects behavioural data into
% a struct that can be loaded for further analyses later, such as plotting 
% or chi square tests. The behav struct has the following organization:

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

clear

%% WHAT DO YOU WANT THIS SCRIPT TO DO

cfg = [];
cfg.output_fd = 'E:\Documents\TmazePaper\data';  % file directory for saving data files
cfg.writeDiary = 0; % keep a text file record of command window history
cfg.writeFiles = 1; % save data
cfg. output_fn = 'behavior'; % what to call the output

%% check that required things exist before continuing 

cfg.rats = {'R042','R044','R050','R064'};% recommended not to change the rats included, do all of them at the same time excepts when testing with the intention of deleting the test data
cfg_temp = []; 
cfg_temp.requireExpKeys = 1;
cfg_temp.ExpKeysFields = {'badTrials'};
cfg_temp.requireMetadata = 1;
cfg_temp.MetadataFields = {'taskvars'};
proceed = checkTmazeReqs(cfg_temp); % make sure we have everything

%% get sessions we can use
if proceed
    
    if cfg.writeDiary % save command window text
        warning off % I don't care, Navi
        cd(cfg.output_fd)
        diary(['behavior','.txt'])
        cd(cfg.output_fd)
        disp(date)
    end
    
    cfg.sessions = {'food','water'};
    originalFolder = pwd;
    
    cfg_temp = []; 
    
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
end

if cfg.writeFiles
    cd(cfg.output_fd)
    save('behavior','behav')
end

if cfg.writeDiary; warning on; diary off; end

cd(originalFolder)