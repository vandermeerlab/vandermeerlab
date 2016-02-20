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
% output struct is called preRate for "precandidate rate" and has the
% following organization:
%
%  precRate.R042.pre = 
%  precRate.R042.taskrest =
%  precRate.R042.post =
%  precRate.R044.pre = 
%  ...
%  preRate.R064.post = 
%
%  The rates are the average rates over all sessions. So, the rates for
%  each session were calculated, and the output is the average rate across
%  6 sessions.

% aacarey Sept 2015

%%

cfg.output_fd = 'E:\Documents\TmazePaper\data';

cfg.writeFiles = 1; % 1 save the data; 0 don't


%% 

iWasHere = pwd;

rats = {'R042','R044','R050','R064'};

all_ratePre = 0;
all_rateTaskRest = 0;
all_ratePost = 0;

% work through each rat
for iRat = 1:length(rats)
    cfg_temp.rats = rats(iRat);
    fd = getTmazeDataPath(cfg_temp);
    
    ratePre = 0;
    rateTaskRest = 0;
    ratePost = 0;
   
    % work through each session folder and do the thing
    for iSession = 1:length(fd)
        cd(fd{iSession})
        LoadExpKeys % for task epoch times
        LoadMetadata % for metadata.taskvars.rest_iv
        evt = load(FindFile('*-candidates.mat')); % the variable this loads is called precandidates too
        % hack to unpack the struct...or whatever
        name = fieldnames(evt);
        precandidates = evt.(name{1});
        
        precPre = restrict(precandidates,ExpKeys.prerecord(1),ExpKeys.prerecord(2)); % prerecord precandidates
        tPre = (ExpKeys.prerecord(2) - ExpKeys.prerecord(1))/60; % length, in minutes, of the prerecord
        ratePre = length(precPre.tstart)/tPre + ratePre; 
        
        precTaskRest = restrict(precandidates,metadata.taskvars.rest_iv.tstart,metadata.taskvars.rest_iv.tend);
        tTaskRest = (sum(metadata.taskvars.rest_iv.tend - metadata.taskvars.rest_iv.tstart))/60;
        rateTaskRest = length(precTaskRest.tstart)/tTaskRest + rateTaskRest;
        
        precPost = restrict(precandidates,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
        tPost = (ExpKeys.postrecord(2) - ExpKeys.postrecord(1))/60;
        ratePost = length(precPost.tstart)/tPost + ratePost;  
    end
    % get avg rate, divide by nSessions
    ratePre = ratePre/length(fd);
    rateTaskRest = rateTaskRest/length(fd);
    ratePost = ratePost/length(fd);
    
    all_ratePre = ratePre + all_ratePre;
    all_rateTaskRest = rateTaskRest + all_rateTaskRest;
    all_ratePost = ratePost + all_ratePost;
    
    % add it to struct (precRate is short for precandidate rate)
    
    precRate.(rats{iRat}).pre = ratePre;
    precRate.(rats{iRat}).taskrest = rateTaskRest;
    precRate.(rats{iRat}).post = ratePost;    
end
% get avg for all rats
precRate.all.pre = all_ratePre/length(rats);
precRate.all.taskrest = all_rateTaskRest/length(rats);
precRate.all.post = all_ratePost/length(rats);


% save that stuff
if cfg.writeFiles
    cd(cfg.output_fd)
    save('SWRMUArate','precRate')
end
disp('End of script run')
cd(iWasHere)