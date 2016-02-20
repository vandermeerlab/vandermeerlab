%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%        Get rate of SWR-associated MUA during session epochs         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Tmaze data
%
% collect number of SWR-associated MUA regions and calculate the rate for
% prerecord, task rest periods, and postrecord. 
% unit is number of precandidates per minute.
%
%  The rates are the average rates over all sessions. So, the rates for
%  each session were calculated, and the output is the average rate across
%  6 sessions.

% aacarey Jan 2015 from poster code

clearvars -except CFG

%%

cfg.output_fd = [pwd,'\visuals'];

cfg.writeFiles = 1; % 1 save the data; 0 don't

cfg.whichCandidates = '-candidates';

cfg.whichITI = 'equalBehavior'; % 'equalBehavior' or 'all'

cfg.output_fn = 'EventRate';

cfg.colormode = 'inventory4'; 

cfg.FontSize = 8;

%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.plot);
end

%% COLLECT DATA

% check spelling
if ~any(strcmp(cfg.whichITI,{'equalBehavior','all'}))
    error('Unrecognized cfg.whichITI option specified')
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~               CALCULATING EVENT RATES                       ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp(' ')

originalFolder = pwd;

% get rat list
rats = TmazeRats;

% categories for restricting event times
categories = {'prerecord',[cfg.whichITI,'ITI'],'task','postrecord'};
restrictionType = {'food','water'};

% initialize variables for all-rat rate
for iCat = 1:length(categories)
    EventRate.all.(categories{iCat}).food = [];
    EventRate.all.(categories{iCat}).water = [];
end

% work through each rat
for iRat = 1:length(rats)
    fprintf('\t*** Working on %s\n',rats{iRat});
    
    cfg_temp.rats = rats(iRat);
    fd = getTmazeDataPath(cfg_temp); % get list of directories for this rat
    
    % initialize variables for current rat
    for iCat = 1:length(categories)
        EventRate.(rats{iRat}).(categories{iCat}).food = 0;
        EventRate.(rats{iRat}).(categories{iCat}).water = 0;
    end
   
    % work through each session folder
    for iSession = 1:length(fd)
        cd(fd{iSession})
        [~,sessionID,~] = fileparts(pwd);
        fprintf('\t\t%s\n',sessionID)
        LoadExpKeys
        LoadMetadata
        evt = loadpop([sessionID,cfg.whichCandidates,'.mat']); evt_original = evt;
        
        %~~~~~~~~~ loop through each event restriction category ~~~~~~~~~~~
        for iCat = 1:length(categories)
            
            %~~~~~~~~~~~    restrict according to category    ~~~~~~~~~~~~~  
            switch categories{iCat}
                case 'all'
                    % keep the whole set
                case 'prerecord'
                    IV = iv(ExpKeys.prerecord(1),ExpKeys.prerecord(2));
                    evt = restrict(evt_original,IV);
                case 'postrecord'
                    IV = iv(ExpKeys.postrecord(1),ExpKeys.postrecord(2));
                    evt = restrict(evt_original,IV);
                case 'task'
                    IV = iv(ExpKeys.task(1),ExpKeys.task(2));
                    evt = restrict(evt_original,IV);
                case 'allITI'
                    IV = metadata.taskvars.rest_iv;
                    evt = restrict(evt_original,IV);
                case 'equalBehaviorITI'
                    IV = GetEqualBehaviorIV([],metadata);
                    evt = restrict(evt_original,IV);
            end % of switch event category
            
            fprintf('\t\t\tnEvents after restricting evt to %s: %d\n',categories{iCat},length(evt.tstart));
            
            if ~isempty(evt.tstart) % happens with equalBehaviorITI
                %~~~~~~~~~ get the length of the interval ~~~~~~~~~~~~~~~~~~~~~
                t =(sum(IV.tend - IV.tstart))/60; % length in minutes
                
                %~~~~~~ Now get session rate and add it to the data collector ~
                EventRate.(rats{iRat}).(categories{iCat}).(ExpKeys.RestrictionType) = length(evt.tstart)/t + EventRate.(rats{iRat}).(categories{iCat}).(ExpKeys.RestrictionType);
            end
            
        end % of categories        
    end % of session
    
    for iCat = 1:length(categories)
        for iRestr = 1:length(restrictionType)
            %~~~~ get avg rate for current rat, divide by nSessions ~~~~~~~~~~~
            EventRate.(rats{iRat}).(categories{iCat}).(restrictionType{iRestr}) = EventRate.(rats{iRat}).(categories{iCat}).(restrictionType{iRestr})/3;
            
            %~~~~ add this rat's rate to the all-rat rate ~~~~~~~~~~~~~~~~~~~~~
            EventRate.all.(categories{iCat}).(restrictionType{iRestr}) = [EventRate.all.(categories{iCat}).(restrictionType{iRestr}) EventRate.(rats{iRat}).(categories{iCat}).(restrictionType{iRestr})];
        end
    end
    
end % of rats

for iCat = 1:length(categories)
    for iRestr = 1:length(restrictionType)
        % get SEM values
        ysemSubj.all.(categories{iCat}).(restrictionType{iRestr}) = nanstd(EventRate.all.(categories{iCat}).(restrictionType{iRestr})/sqrt(length(rats)));
        %~~~~~~~~~~~~ get avg for all rats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        EventRate.all.(categories{iCat}).(restrictionType{iRestr}) = nansum(EventRate.all.(categories{iCat}).(restrictionType{iRestr}))/length(rats);
    end
end

%% PLOTTING SECTION

colors = TmazeColors('inventory4');

epochs = categories; % because I call them epochs in other plotting scripts
rats = ['all',TmazeRats];

location = [1 2.5]; % where to place the bars
xlims = [0 location(end)+1];

position_all = [1 2 9 10; 3 4 11 12; 5 6 13 14; 7 8 15 16];
start = [17 18 25 26];
jhgjhg = [0 2 4 6];

figure; hold on

for iEpoch = 1:length(epochs)
    d = [];
    col = {};
    legtxt = {};
    for iRat = 1:length(rats)
        
        switch rats{iRat}
            case 'all'
                position = position_all(iEpoch,:);
                txt = 'all rats';
            otherwise
                position = start(iRat-1) + jhgjhg(iEpoch); % tracks which spot the subplot will occupy
                txt = rats{iRat};
        end
              
        % get colors into correct format
        col = {colors.(rats{iRat}).f colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
        
        % get data into format for plotting
        d = [EventRate.(rats{iRat}).(epochs{iEpoch}).food EventRate.(rats{iRat}).(epochs{iEpoch}).water];
     
    % for ylims assignment
    extreme(iEpoch,iRat) = max(d);
    
    h.(rats{iRat})(iEpoch) = subplot(4,8,position);
        for iBar = 1:length(d)
            bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
            hold on
        end
    
         if isequal(position,[1 2 9 10]) || position(1) == 17 || position(1) == 25
             if iRat ==1 
            ylab = 'Average rate of events (n/min)';
             else
               ylab = {'Average rate of ','events (n/min)'};  
             end
             ylabel(ylab,'FontSize',cfg.FontSize)
        else
            set(gca,'YTick',[]);
         end
     
         if iRat == 1
             switch epochs{iEpoch}
                 case 'prerecord'
                     title('PRERECORD','FontSize',cfg.FontSize)
                 case 'postrecord'
                     title('POSTRECORD','FontSize',cfg.FontSize)
                 case 'task'
                     title('TASK','FontSize',cfg.FontSize)
                 case 'allITI'
                     title('INTERTRIAL','FontSize',cfg.FontSize)
                 case 'equalBehaviorITI'
                     title('EQUAL BEHAVIOR ITI','FontSize',cfg.FontSize)
             end
             
             errorbar(location,d,[ysemSubj.all.(epochs{iEpoch}).food ysemSubj.all.(epochs{iEpoch}).water],'Color','k','LineStyle','none');
             text(0.62,0.9,txt,'Units','normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','right',...
                'FontSize',cfg.FontSize)
         else
             text(0.7,0.85,txt,'Units','normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','right',...
                'FontSize',cfg.FontSize)
         end
             
    set(gca,'XLim',xlims); box off
    end
end % of epochs

for iRat = 1:length(rats)
    set(h.(rats{iRat}),'XTick',location,'XTickLabel',{'food' 'water'},'FontSize',cfg.FontSize,'LineWidth',1,'XLim',xlims,'Layer','top');
    
    for iEpoch = 1:length(epochs)
        
            set(h.(rats{iRat})(iEpoch), 'YLim', [0 50])
     
    end % of epochs
end % of rats

hold off

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

% save it
if cfg.writeFiles
    cd(cfg.output_fd)
    print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
end

disp(' ')
disp('~~~~~~~ END OF EVENT RATE RUN ~~~~~~~')
cd(originalFolder)