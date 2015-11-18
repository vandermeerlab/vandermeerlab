%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% Get data that will be used for Co-Occurrence and Sequence Analysis  %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Tmaze data
% Collect data relevant to co-occurrence and sequence analysis into a
% single struct. 
% The collected input data goes into a struct called "inputData": this is any
% data that is used to get results, and contains data from ALL FOUR RATS at
% the same time, so we know that the input data was generated using the same
% parameters.

% aacarey Sept 2015

clear 

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

cfg.writeFiles = 1; % if 1, saves inputData as outlined below; if 0, does not save inputData

% what to call the output and where to put it?
cfg.output_fn = 'inputData';
cfg.output_fd = pwd; % 'E:\Documents\TmazePaper\data';

% Which cells to load for LoadSpikes?
cfg.load_questionable_cells = 1; % 1 load 5's, 0 don't

cfg.maxspeed = 3.5; % in cm/s, upper speed limit for estimating tuning curves

%%
please = []; please.checkall = 1;
if 0 %~checkTmazeReqs(please)
    return;
end

tic

iWasHere = pwd; % remember where I started before cding all over the place

% don't change this...keep all rats in the inputData. Can choose in
% later scripts to take one rat or another
rats = {'R042','R044','R050','R064'}; % 'R042','R044','R050','R064'
arms = {'L','R'}; % shortforms for left and right arms, used for storing corresponding data

for iRat = 1:length(rats)
    please.rats = rats(iRat);
    fd = getTmazeDataPath(please); % get list of directories for session folders
    
    disp(' '); disp('~~~~~~~~~~~~'); disp(['~~  ',rats{iRat},'  ~~']); disp('~~~~~~~~~~~~')
    
    % work through each session
    for iSession = 1:length(fd)
        cd(fd{iSession})
        [~,sessionID,~] = fileparts(pwd);
        disp(' '); disp(['Working on ',sessionID])
        
        % load the data
        
        LoadMetadata; LoadExpKeys
        
        please = [];
        please.load_questionable_cells = cfg.load_questionable_cells;
        S = LoadSpikes(please);
        
        please = [];
        please.convFact = ExpKeys.convFact; % load position in units of cm
        pos = LoadPos(please);
        
        % save some info about this session
        sessionData.sessionID = sessionID;
        sessionData.restrictionType = ExpKeys.RestrictionType;
        
        % deal with trial intervals
        if strcmp(rats(iRat),'R042')
            metadata = TrimTrialTimes([],metadata); % R042 only!!
        end
        
        [iv_L,iv_R] = GetMatchedTrials([],metadata,ExpKeys); % get equal number of L and R trials
        
        % put trial intervals into the data collector
        sessionData.L.iv = iv_L;
        sessionData.R.iv = iv_R;
        
        % put coords into the data collector
        sessionData.L.coord = metadata.coord.coordL_cm;
        sessionData.R.coord = metadata.coord.coordR_cm;
        please = [];
        please.binSize = 3;
        please.run_dist = ExpKeys.pathlength;
        for iArm = 1:length(arms)
            please = [];
            please.binSize = 3;
            please.run_dist = ExpKeys.pathlength;
            % my coord is already in cm, but this will put it into new bin sizes
            sessionData.(arms{iArm}).coord = StandardizeCoord(please,sessionData.(arms{iArm}).coord);
        end
        
        chp = tsd(0,metadata.coord.chp_cm,{'x','y'}); % make choice point useable by cobebase functions
        
        % Get intervals where the rat is running quickly. The point is
        % to exclude times when he is going so slow that he's likely
        % making SWRs (participating cells are spiking non-locally,
        % which affects tuning curve estimates).
        
        spd = getLinSpd([],pos);
        
        please = []; please.method = 'raw'; please.threshold = cfg.maxspeed; % in cm/sec
        iv_run = TSDtoIV(please,spd); % intervals with speed above 3.5 cm/s
        
        
        % restrict spike and position data
        
        for iArm = 1:length(arms)
            % first, keep full (unrestricted) S
            sessionData.(arms{iArm}).S_full = S;
            
            % restrict to L or R intervals only
            sessionData.(arms{iArm}).S = restrict(S,sessionData.(arms{iArm}).iv);
            % on the first pass, this ^^ reads as sessionData.L.S = restrict(S,sessionData.L.iv)
            sessionData.(arms{iArm}).pos = restrict(pos,sessionData.(arms{iArm}).iv);
            
            % restrict to fast intervals only
            sessionData.(arms{iArm}).pos = restrict(pos,iv_run);
            sessionData.(arms{iArm}).S = restrict(S,iv_run);
        end
        
        % linearize position data and choice point
        for iArm = 1:length(arms)
            please = []; please.Coord = sessionData.(arms{iArm}).coord;
            sessionData.(arms{iArm}).linpos = LinearizePos(please,sessionData.(arms{iArm}).pos);
            sessionData.(arms{iArm}).chp = LinearizePos(please,chp); % the exact chp for L or R differs depending on the coord
            
        end
        
        % get the tuning curves and fields for L or R trajectory
        for iArm = 1:length(arms)
            please = [];
            please.binSize = 1; % if 1, keep same 3 cm bins as were used in StandardizeCoord
            sessionData.(arms{iArm}).PF = MakeTC(please,sessionData.(arms{iArm}).S,sessionData.(arms{iArm}).linpos);
        end
        
        % get cells that are unique to a single arm (rather than whole
        % trajectory)
        % remove cells with fields before choice point
        for iArm = 1:length(arms)
            fields.(arms{iArm}) = [];
            fields.(arms{iArm}) = sessionData.(arms{iArm}).PF.field_template_idx(sessionData.(arms{iArm}).PF.field_loc > sessionData.(arms{iArm}).chp.data(1));
        end
        % find unique cells (L or R only)
        [~,removeL,removeR] = intersect(fields.L,fields.R);
        fields.L(removeL) = []; fields.R(removeR) = [];
        
        % make a new field for unique cells and order them
        for iArm = 1:length(arms)
            % unique spike trains
            sessionData.(arms{iArm}).Sunique = OrderSelectS([],sessionData.(arms{iArm}).S,(fields.(arms{iArm})));
            sessionData.(arms{iArm}).Sunique_full = OrderSelectS([],sessionData.(arms{iArm}).S_full,(fields.(arms{iArm})));
        end
        
        % order the non-unique cells too
        for iArm = 1:length(arms)
            sessionData.(arms{iArm}).S = OrderSelectS([],sessionData.(arms{iArm}).S,(sessionData.(arms{iArm}).PF.field_template_idx));
            sessionData.(arms{iArm}).S_full = OrderSelectS([],sessionData.(arms{iArm}).S_full,(sessionData.(arms{iArm}).PF.field_template_idx));
            %sessionData.(arms{iArm}).S.t = S.t(sessionData.(arms{iArm}).PF.template_idx);
            %sessionData.(arms{iArm}).S.label = S.label(sessionData.(arms{iArm}).PF.template_idx);
            %sessionData.(arms{iArm}).S.usr.data = S.usr.data(sessionData.(arms{iArm}).PF.template_idx);
        end
        
        % make useful tuning curves that are ordered
        for iArm = 1:length(arms)
            please = [];
            sessionData.(arms{iArm}).tc_ordered = sessionData.(arms{iArm}).PF.tc(sessionData.(arms{iArm}).PF.field_template_idx,:);
            sessionData.(arms{iArm}).tc_unique_ordered = sessionData.(arms{iArm}).PF.tc(fields.(arms{iArm}),:);
        end
        
        all_sessionData.(rats{iRat})(iSession) = sessionData;
    end
    
    
    
    %all_sessionData.(rats{iRat}).(['s',num2str(iSession)]) = sessionData;
end

% make output struct
inputData = all_sessionData;
inputData.date = datestr(now);

if cfg.writeFiles
    disp(['Writing inputData to .mat file in ',cfg.output_fd,'...'])
    cd(cfg.output_fd)
    save(cfg.output_fn,'inputData')
else
    disp('WARNING: You have selected not to save the inputData')
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~            End of inputData run                ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
toc
cd(iWasHere)
