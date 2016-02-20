%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%      Get Number of Events within Different Restrict Categories      %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script works through session data folders and counts the number of
% events that occur within restrict categories that we employ throughout
% other analysis steps. These restrict categories are not to be confused
% with restriction type. Rather, these are intervals that we define before
% analysis that allow us to group and compare the content of spike or LFP 
% data during detected events. 
% Examples would be restricting our detected SWR events to the prerecord and
% postrecord phases of the experiment so that we can compare the content
% of SWR-associated neural activity before and after task performance.
% This script was written to address the question of how many events are
% represented in each restrict category.
%
% The output variable is called nEvt and has the following fields:
%
% nEvt.R042 ...
% nEvt.R044 ...
% nEvt.R050 ...
% nEvt.R064 ...
%
% Within each rat ID field are additional fields corresponding to the
% restriction type ('food' and 'water') as well as a 'total' field for all
% events detected across 6 session within a certain category, regardless of
% restriction type. 
% The food and water fields contain data accumulated from 3 sessions with 
% corresponding restriction types.
%
% Typically, food days are sessions 1, 3, 5 and water days are 2, 4, 6.
%
% nEvt.R042.food ...
% nEvt.R042.water ...
% nEvt.R042.total ... food and water counts combined
%
% Each restriction field contains 6 additional fields corresponding to
% restrict categories (intervals that we can restrict or limit our detected
% events to):
%
% nEvt.R042.food.all = N; % the total number of events detected on food days
% nEvt.R042.food.prerecord = n; % the number of prerecord events on food days
% nEvt.R042.food.task = n; % task events on food days
% nEvt.R042.food.postrecord = n; % postrecord events on food days
% nEvt.R042.food.allITI = n; % intertrial events on food days (intertrial
%                         periods are those during which the rat was sitting 
%                         on a waiting platform rather than performing the
%                         task on the Tmaze.
% nEvt.R042.food.equalBehaviorITI = n; 
%
% The nEvt variable is output for potential use in plotting or table-generating
% scripts down the line. Otherwise, the text file or contents of the nEvt
% variable can serve as a direct reference for event counts reported in the
% paper.

% aacarey Dec 1 2015

clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% Which candidates file?
cfg.whichCandidates = '-candidates'; % the suffix of the candidates file you want to load

% Would you like to save the output?
cfg.writeFiles = 1; % If 1, save the output variable as well as a record of
% command window history. If 0, don't.

% Where to save the output?
cfg.output_fd = [pwd,'\data'];

% What to call the output?
cfg.output_fn = 'nEvents';

% Which rats to use? (best to leave all here)
cfg.rats = TmazeRats;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                       Main body of script                           %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic % start timer

% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.plot);
end


iWasHere = pwd; % remember original directory

if cfg.writeFiles % start the text record of command window history
    cd(cfg.output_fd)
    diary([cfg.output_fn,'_text.txt'])
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                  COLLECTING nEVENT DATA                     ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

%~~~~~~~~~~~~~~~~~~~~~~~~~ for the record ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf('\nDATE of script run: %s\n',datestr(now));
fprintf('COMPUTER of data origin: %s\n',getenv('COMPUTERNAME'));

%~~~~~~~~~~~~~~~~~~ declare event restriction categories ~~~~~~~~~~~~~~~~~~
categories = {'all','prerecord','task','postrecord','allITI','equalBehaviorITI'};

%~~~~~~~~~~~~~~~~~~~~~~ iterate through each rat ID ~~~~~~~~~~~~~~~~~~~~~~~
for iRat = 1:length(cfg.rats)
    disp(' '); disp('~~~~~~~~~~~~~~~~~~~~~~~~'); disp(['~~        ',cfg.rats{iRat},'        ~~']); disp('~~~~~~~~~~~~~~~~~~~~~~~~')
    
    % initialize some empty substructs
    nEvt.(cfg.rats{iRat}).food = [];
    nEvt.(cfg.rats{iRat}).water = [];
    
    % ask for a list of directories for the current rat's session data
    cfg_temp.rats = cfg.rats(iRat);
    fd = getTmazeDataPath(cfg_temp); % session directories for current rat        
    
    %~~~~~~~~~~~~~~~~~~~~ work through each session ~~~~~~~~~~~~~~~~~~~~~~~
    for iSession = 1:length(fd)
        
        cd(fd{iSession})
        [~,sessionID,~] = fileparts(pwd);
        disp(' '); disp([' *** Working on ',sessionID,' ***'])
        
        %~~~~~~~~~~~~~~~~~~~~~~~ Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ExpKeys = []; LoadExpKeys; disp(['     Restriction type: ',ExpKeys.RestrictionType])
        evt = []; evt = loadpop([sessionID,cfg.whichCandidates,'.mat']); evt_original = evt;
        metadata = []; LoadMetadata
        
        % tell me how many total events were detected
        fprintf('      nEvents loaded: %d\n',length(evt.tstart));
        
        %~~~~~~~~~ loop through each event restriction category ~~~~~~~~~~~
        for iCat = 1:length(categories)
            
            %~~~~~~~~~~~    restrict according to category    ~~~~~~~~~~~~~  
            switch categories{iCat}
                case 'all'
                    % keep the whole set
                case 'prerecord'
                    evt = restrict(evt_original,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
                case 'postrecord'
                    evt = restrict(evt_original,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
                case 'task'
                    evt = restrict(evt_original,ExpKeys.task(1),ExpKeys.task(2));
                case 'allITI'
                    evt = restrict(evt_original,metadata.taskvars.rest_iv.tstart(1),metadata.taskvars.rest_iv.tend(end));
                case 'equalBehaviorITI'
                    eq_iv = GetEqualBehaviorIV([],metadata);
                    evt = restrict(evt_original,eq_iv);
            end % of switch event category
            
            fprintf('       nEvents after restricting evt to %s: %d\n',categories{iCat},length(evt.tstart));
            
            %~~~~~~~~~~~~~~~~~~ sum the numbers ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            % check if the value has been initialized, if not, then start
            % sum at 0
            if ~isstruct(nEvt.(cfg.rats{iRat}).(ExpKeys.RestrictionType)) || ~isfield(nEvt.(cfg.rats{iRat}).(ExpKeys.RestrictionType),categories{iCat})
                nEvt.(cfg.rats{iRat}).(ExpKeys.RestrictionType).(categories{iCat}) = 0;
            end
            % add event count to existing count
            nEvt.(cfg.rats{iRat}).(ExpKeys.RestrictionType).(categories{iCat}) = nEvt.(cfg.rats{iRat}).(ExpKeys.RestrictionType).(categories{iCat}) + length(evt.tstart);
                       
        end % of categories loop 
        
        %~~~~~~~~ get total counts for both restr types together ~~~~~~~~~~
        
        if iSession == length(fd) % sum everything
            for iCat = 1:length(categories)
                nEvt.(cfg.rats{iRat}).total.(categories{iCat}) = nEvt.(cfg.rats{iRat}).food.(categories{iCat}) + nEvt.(cfg.rats{iRat}).water.(categories{iCat});
            end
            
            % tell me the counts
            totalTypes = {'food','water','total'};
            for iType = 1:length(totalTypes)
                disp(' '); disp([' *** Displaying ',totalTypes{iType},' counts for ',cfg.rats{iRat},' ***'])
                for iCat = 1:length(categories)
                    fprintf('       %s: %d \n',categories{iCat}, nEvt.(cfg.rats{iRat}).(totalTypes{iType}).(categories{iCat}))
                end
            end
        end
    end % of session loop
end % of rat loop

%%
disp(' ')
if cfg.writeFiles % save the output
    disp(['Writing nEvents data to .mat file in ',cfg.output_fd,'...'])
    cd(cfg.output_fd)
    save(cfg.output_fn,'nEvt')
else
    disp('WARNING: You have selected not to save the output data')
end

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                End of nEvents run                ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

toc % stop timer

if cfg.writeFiles; diary off; end
cd(iWasHere) % return to original directory