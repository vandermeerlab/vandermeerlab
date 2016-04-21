%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% Get data that will be used for Co-Occurrence and Sequence Analysis  %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Tmaze paper
% Collect data relevant to co-occurrence and sequence analysis into a
% single struct. 
% The collected input data goes into a struct called "inputData": this is any
% data that is used to get results, and contains data from ALL FOUR RATS at
% the same time, so we know that the input data was generated using the same
% parameters.
%
% inputData and its corresponding text file (record of command window
% history) is saved to the current directory unless otherwise specified.
%
% *** some script behaviors are configurable, see the section titled "WHAT
% DO YOU WANT THIS SCRIPT TO DO?" and read the descriptions.
%
% ORGANIZATION OF INPUTDATA
% The inputData collector is a 1x1 struct with 6 fields. The first field is
% a header with information about the particular script run that produced
% the inputData
%
% inputData.hdr.date = datestr(now); % the date and time when the script was
%             run, to serve as a type of identity
% inputData.hdr.version = version; % the version of MATLAB that ran the
%             script
% inputData.hdr.mfilename = mfilename; % the name of this script
% inputData.hdr.host = 'CALLISTO'; % the name of the machine that ran the
%             script
%
% inputData also saves the script's config in inputData.cfg = cfg;
%
% % additional inputData fields correspond to the data coming from each rat:
% inputData.R042...
% inputData.R044...
% inputData.R050...
% inputData.R064...
%
% % within each rat ID field is a 1x6 struct with 4 fields (see below).
% Each of the 6 structs corresponds to a recording day. inputData.R042(1)
% is recording day 1 and so on to R042(6) which is recording day 6.
%
% inputData.R042(1).sessionID = 'R042-2013-08-16'; % the date of data
%                             collection / the id of the source folder data
% inputData.R042(1).restrictionType = 'food';
% inputData.R042(1).L... ; the data corresponding to the L arm
% inputData.R042(1).R... ; the data corresponding to the R arm
% 
% % Within each L and R field is a 1x1 struct with 12 fields (at the time
% this was written):
%
% inputdata.R042(1).L.iv =...; all the trial intervals that are used to
%                    estimate tuning curves, the output from
%                    GetMatchedTrials
%
%                ...L.coord = ...; the standardized coord corresponding to
%                    the left trajectory
%
%                ...L.S_restr = ...; the S restricted to left trials when
%                    the rat was running quickly enough that we expect him
%                    not to be making SWRs (for estimating tuning curves).
%
%                ...L.pos_restr = ...; as with S above, but for position
%                    data
%
%                ...L.linpos = ...; linearized pos_restr, the output from
%                   LinearizePos
%
%                ...L.chp = ...; the location of the choice point in
%                    linearized position data
%
%                ...L.PF = ...; place field information and tuning curves,
%                    the output from MakeTC
%
%                ...L.S_traj = ...; all of the units that have place fields
%                    along the left trajectory, ordered according to place field
%                    location from the start to the end of the left arm
%
%                ...L.S_arm = ...; all of the units that have place fields on
%                    the left arm, ordered (these will be used for sequence
%                    analysis)
%
%                ...L.S_arm_unique = ...; all of the units that have place
%                    fields on the left arm only, ordered; no fields on the right
%                    arm (these will be used for co-occurrence analysis)
%
%                 ...L.tc_traj = ...; ordered tuning curves for the
%                    units in S_traj
%
%                 ...L.tc_arm = ...; ordered tuning curves for the
%                    units in S_arm
%
%                 ...L.tc_arm_unique = ...; ordered tuning curves for
%                    the units in S_arm_unique
%
% aacarey Sept 2015, Nov 2015

clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% save the data?
cfg.writeFiles = 1; % if 1, saves inputData and text file as described below; if 0, does not save 

% what to call the output and where to put it?
cfg.output_fn = 'inputData'; % output filename. will also save a text file with the name (for example) inputData_text
cfg.output_fd = [pwd,'\data']; % the output file directory, ex 'E:\Documents\TmazePaper\data';

% Which cells to load for LoadSpikes?
cfg.load_questionable_cells = 1; % 1 load 5's, 0 don't

% Which position data units to work with? (For testing)
cfg.posunits = 'px'; % 'cm' for centimeters or 'px' for pixels. Note that either way, the data gets translated into cm

switch cfg.posunits
    case 'px'
        cfg.maxspeed = 10; % in pixels per second
    case 'cm'
        cfg.maxspeed = 3.5; % in cm/s, upper speed limit for estimating tuning curves
end

% remove bad trials? (these are trials where the rat attempted to reverse
% directions and/or was interfered with by the researcher)
cfg.rmbadtrl = 1; % 1 yes, 0 no

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                          DO THE THING                               %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *** None of the following code should be altered between runs. If additional
% configurations are needed, add them to the section above.

% This script uses dynamic fieldnames. What this means that, for example, I
% can iterate through a cell array of desired "fieldnames" and assign them
% to a struct through different passes of the same loop. In this way I can 
% be sure the same config and function is being applied to L and R data. 
% Here's a brief example:

% fieldnames = {'L','R'};
% 
% collect.question = 'How do I win this game of chess?';
% 
% for iFieldName = 1:length(fieldnames)
%     collect.(fieldnames{iFieldName}).answer = ['Just take a ',fieldnames{iFieldName},'ook, it''s easy.'];
% end

% now if i look inside of the variable collect, I have the following
% organization:

% collect.L.answer = 'Just take a Look, it's easy.';
% collect.R.answer = 'Just take a Rook, it's easy.';

% So basically, since both L and R are going through the same loop, I can
% know that the only difference between the output strings is that fact
% that I use L or R. If I kept them separate (did them outside of a loop)
% then I have a greater risk of accidentally using a different parameter
% somewhere. It's a *lot* harder to read and interpret, but it's makes the
% analysis more "secure".

% Note that dynamic fieldnames are not the only way to do this, but I like
% the easy way of being able to see the fieldnames and opening the data
% collector in the Variables tab in MATLAB.
% In the case of recording day 1 data and so on (the 1x6 structs in each
% rat ID field of inputData), these are done using dynamically changing 1xn
% structs through each iteration.

%% proceed with other things

iWasHere = pwd; % remember where I started before cding all over the place

%~~~~~~~~~~~ handle saving a record of command window text ~~~~~~~~~~~~~~~~
if cfg.writeFiles
    cd(cfg.output_fd)
    diary([cfg.output_fn,'_text.txt'])
    cd(iWasHere)
end

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                    Collecting inputData                     ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
fprintf('\nDATE of script run: %s\n',datestr(now));

%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.in);
    % and we already checked requisites
else
    %~~~~~~~~~~~~~~ verify that all requisite data exists ~~~~~~~~~~~~~~~~~~~~~
    please = []; please.checkall = 1;
    if ~checkTmazeReqs(please); return; end
end

tic

%~~~~~~~~~~~~~ save some info about this particular run ~~~~~~~~~~~~~~~~~~~
inputData.hdr.date = datestr(now); % the datestr can serve as a data "ID" because no two runs will have exactly the same date and time
inputData.hdr.version = version; % the version of matlab than ran the script
inputData.hdr.mfilename = mfilename; % the name of this script (can help point to the version, in some cases)
disp (' '); disp(['MFILENAME: ',inputData.hdr.mfilename])
inputData.hdr.host = getenv('COMPUTERNAME'); % the computer's name that ran the script (if we know the date and filename, 
% and GitHub history for this machine, then we can know even more about the version of the script that produced the inputData)
disp(' '); disp(['COMPUTERNAME: ',inputData.hdr.host])

%~~~~~~~~~~~~~~~~~~~~~ set the dynamic field names ~~~~~~~~~~~~~~~~~~~~~~~~
rats = {'R042','R044','R050','R064'}; % do not change from 'R042','R044','R050','R064'
arms = {'L','R'}; % shortforms for left and right arms, used for storing corresponding data

nArms = length(arms);

%~~~~~~~~~~~~~~~~~~~~~~ iterate through each rat ID ~~~~~~~~~~~~~~~~~~~~~~~
% this is the loop that generates, for example, inputData.R042
for iRat = 1:length(rats)
    please = [];
    please.rats = rats(iRat);
    fd = getTmazeDataPath(please); % get list of directories for session folders
    
    disp(' '); disp('~~~~~~~~~~~~~~~~'); disp(['~~    ',rats{iRat},'    ~~']); disp('~~~~~~~~~~~~~~~~')
    
    %~~~~~~~~~~~~~~~~~~~~ work through each session ~~~~~~~~~~~~~~~~~~~~~~~
    % this is the loop that generates the 1x6 struct, for example
    % inputData.R042(1) and inputData.R042(2) and so on
    for iSession = 1:length(fd)
        cd(fd{iSession})
        [~,sessionID,~] = fileparts(pwd);
        disp(' '); disp(['*** Working on ',sessionID,' ***'])
        
        % ~~~~~~~~~~~~~~~~~~~~ LOAD DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % load the metadata and keys
        metadata = []; ExpKeys = []; % initialize and/or clear existing
        LoadMetadata; LoadExpKeys
        
        disp(['Restriction type: ',ExpKeys.RestrictionType]);
        
        % load spiking data
        please = [];
        please.load_questionable_cells = cfg.load_questionable_cells;
        please.getTTnumbers = 1;
        please.getRatings = 1;
        S = LoadSpikes(please);
        
        % load LFP data
        please = [];
        please.fc = ExpKeys.goodSWR(1);
        CSC = LoadCSC(please);
        
        % add a usr firing rate field in S
        please = [];
        S = AddRateTS(please,S,CSC.tvec);
        
        % add a usr field that allows cell# tracking when cells are selected and reordered
        S.usr.cell_num = 1:length(S.t);
        
        % load position data
        please = [];
        if strcmp(cfg.posunits,'cm'); please.convFact = ExpKeys.convFact; end % load position in units of cm, if requested
        cfg.removeZeros = 1; cfg.tsflag = 'sec';
        pos = LoadPos(please);
        
        %~~~~~~~~~~~~~ save some info about this session ~~~~~~~~~~~~~~~~~~
        sessionData.sessionID = sessionID;
        sessionData.restrictionType = ExpKeys.RestrictionType;
        
        %~~~~~~~~~~~~~~~~~ deal with trial intervals ~~~~~~~~~~~~~~~~~~~~~~
        
        % fix R042's trial times to match the chewing restriction times
        if strcmp(rats(iRat),'R042')
            disp('Trimming trials for R042')
            metadata = TrimTrialTimes([],metadata); % R042 only!!
        end
        
        % get the same number of L and R trials
        switch cfg.rmbadtrl
            case 0
                [iv_L,iv_R] = GetMatchedTrials_old([],metadata); % get equal number of L and R trials
            case 1
                [iv_L,iv_R] = GetMatchedTrials([],metadata,ExpKeys); % get equal number of L and R trials
            otherwise
                error('Unrecognized cfg.rmbadtrl')
        end
        
        %~~~~~~~~~ put trial intervals into the data collector ~~~~~~~~~~~~
        sessionData.L.iv = iv_L;
        sessionData.R.iv = iv_R;
        
        %~~~~~~~~~~~~~ put coords into the data collector ~~~~~~~~~~~~~~~~~
        % get correct coord, depending on whether user wants to work in
        % units of pixels or centimeters
        switch cfg.posunits
            case 'px'
                sessionData.L.coord = metadata.coord.coordL;
                sessionData.R.coord = metadata.coord.coordR;
            case 'cm'
                sessionData.L.coord = metadata.coord.coordL_cm;
                sessionData.R.coord = metadata.coord.coordR_cm;     
        end
        
        % resample coords to have same number of pos units in each bin
        please = [];
        please.binsize = 3;
        please.run_dist = ExpKeys.pathlength;
        for iArm = 1:nArms
            % my coord is already in cm, but this will put it into new bin sizes
            sessionData.(arms{iArm}).coord = StandardizeCoord(please,sessionData.(arms{iArm}).coord);
        end
        
        % get location of choice point
        switch cfg.posunits
            case 'px'
                chp = tsd(0,metadata.coord.chp,{'x','y'}); % make choice point useable by cobebase functions
            case 'cm'
                chp = tsd(0,metadata.coord.chp_cm,{'x','y'}); % make choice point useable by cobebase functions
        end
        
        %~~~~~~~~~~~~~ restrict spike and position data ~~~~~~~~~~~~~~~~~~
        
        % Get intervals where the rat is running quickly. The point is
        % to exclude times when he is going so slowly that he's likely
        % making SWRs (participating cells are spiking non-locally,
        % which affects tuning curve estimates for local firing). We also
        % restrict the data to left and right trials.
        
        % linear speed
        spd = getLinSpd([],pos);
        
        % find intervals where the linear speed is above threshold
        please = []; please.method = 'raw'; please.threshold = cfg.maxspeed; % in cm/sec
        iv_run = TSDtoIV(please,spd); % intervals with speed above thresh
                
        % restrict spike and position data to trial intervals and running
        % intervals
        disp('~~ Restricting spike and position data to trial intervals and running intervals ~~')
        for iArm = 1:nArms
           
            % restrict to L or R intervals only
            sessionData.(arms{iArm}).S_restr = restrict(S,sessionData.(arms{iArm}).iv);
            % on the first pass, this ^^ reads as sessionData.L.S_restr = restrict(S,sessionData.L.iv)
            sessionData.(arms{iArm}).pos_restr = restrict(pos,sessionData.(arms{iArm}).iv);

            % restrict to fast intervals only
            sessionData.(arms{iArm}).S_restr = restrict(sessionData.(arms{iArm}).S_restr,iv_run);
            sessionData.(arms{iArm}).pos_restr = restrict(sessionData.(arms{iArm}).pos_restr,iv_run); % this line is the main difference between inputData and Generate_CoOccur
        end
        
        %~~~~~~~~~~~ linearize position data and choice point ~~~~~~~~~~~~~
        for iArm = 1:nArms
            please = []; please.Coord = sessionData.(arms{iArm}).coord;
            sessionData.(arms{iArm}).linpos = LinearizePos(please,sessionData.(arms{iArm}).pos_restr);
            sessionData.(arms{iArm}).chp = LinearizePos(please,chp); % the exact chp for L or R differs depending on the coord            
        end
        
        %~~~~~ get the tuning curves and fields for L or R trajectory ~~~~~
        disp('~~ Finding place cells ~~')
        for iArm = 1:nArms
            disp(['  ',arms{iArm},' arm...']); fprintf('\b')
            please = [];
            please.binSize = 1; % if 1, keep same 3 cm bins as were used in StandardizeCoord
            sessionData.(arms{iArm}).PF = MakeTC(please,sessionData.(arms{iArm}).S_restr,sessionData.(arms{iArm}).linpos);
        end
        
        %~~~~~~~~~~~~~~~ order full S for whole trajectory ~~~~~~~~~~~~~~~~
        disp('~~ Ordering all place cells ~~')
        for iArm = 1:nArms
            disp(['  ',arms{iArm},' arm...']); fprintf('\b')
            please = []; please.verbose = 1;
            sessionData.(arms{iArm}).S_traj = SelectTS(please,S,sessionData.(arms{iArm}).PF.template_idx);
        end
        
        %~~~~~~ get cells that were active on only one trajectory ~~~~~~~~~
        disp('~~ Ordering and selecting unique trajectory place cells (fire only for a L or R run) ~~')
        
        % find cell IDs that are common to both L and R spiketrains
        [~,remove.L,remove.R] = intersect(sessionData.L.S_traj.usr.cell_num,sessionData.R.S_traj.usr.cell_num);
        
        % remove the common cells from each group
        for iArm = 1:nArms
            keep.(arms{iArm}) = sessionData.(arms{iArm}).S_traj.usr.cell_num;
            keep.(arms{iArm})(remove.(arms{iArm})) = [];
            disp(['  ',arms{iArm},' arm...']); fprintf('\b')
            please = []; please.verbose = 1;
            sessionData.(arms{iArm}).S_traj_unique = SelectTS(please,S,keep.(arms{iArm}));
        end
        
        %~~~~~~~~~~ find cells with fields after choice point ~~~~~~~~~~~~~
        for iArm = 1:nArms
            fields.(arms{iArm}) = [];
            fields.(arms{iArm}) = sessionData.(arms{iArm}).PF.field_template_idx(sessionData.(arms{iArm}).PF.field_loc > sessionData.(arms{iArm}).chp.data(1)); % oh my
            uniquefields.(arms{iArm}) = []; % initialize as empty for later step
        end
        
        %~~~~~~~~~~~~~~ get cells with fields on the arms ~~~~~~~~~~~~~~~~~
        disp('~~ Ordering and selecting non-unique place cells (may fire on one or both arms) ~~')
        for iArm = 1:nArms
            disp(['  ',arms{iArm},' arm...']); fprintf('\b')
            please = []; please.verbose = 1;
            sessionData.(arms{iArm}).S_arm = SelectTS(please,S,fields.(arms{iArm}));
        end
        
        %~~~~~~~~~~~~~~~ find unique cells (L or R only) ~~~~~~~~~~~~~~~~~~
        % find fields that are common to both groups
        [~,remove.L,remove.R] = intersect(fields.L,fields.R);
        uniquefields = fields;
        % remove the cells with common fields
        for iArm = 1:nArms
            uniquefields.(arms{iArm})(remove.(arms{iArm})) = [];
            
            % do this to maintain order (unique sorts ascending, but we want
            % same order as what went it)
            [~,idx.(arms{iArm})] = unique(uniquefields.(arms{iArm}));
            uniquefields.(arms{iArm}) = sort(idx.(arms{iArm}));
        end
        
        %~~~~~~~~~~~~~~~ get unique cells and order them ~~~~~~~~~~~~~~~~~~
        disp('~~ Ordering and selecting unique place cells (fire on one arm only) ~~')
        for iArm = 1:nArms
            disp(['  ',arms{iArm},' arm...']); fprintf('\b')
            please = []; please.verbose = 1;
            sessionData.(arms{iArm}).S_arm_unique = SelectTS(please,S,(uniquefields.(arms{iArm})));
        end
               
        %~~~~~~~~~~~~~~~~~~~~ order the tuning curves ~~~~~~~~~~~~~~~~~~~~~
        for iArm = 1:nArms
            disp(['Ordering ',arms{iArm},' tuning curves'])
            please = [];
            sessionData.(arms{iArm}).tc_traj = sessionData.(arms{iArm}).PF.tc(sessionData.(arms{iArm}).PF.template_idx,:);
            sessionData.(arms{iArm}).tc_traj_unique = sessionData.(arms{iArm}).PF.tc(keep.(arms{iArm}),:);
            sessionData.(arms{iArm}).tc_arm = sessionData.(arms{iArm}).PF.tc(fields.(arms{iArm}),:);
            sessionData.(arms{iArm}).tc_arm_unique = sessionData.(arms{iArm}).PF.tc(uniquefields.(arms{iArm}),:);
        end
        
        inputData.(rats{iRat})(iSession) = sessionData;
        % inputData.R042(1) = sessionData;
    end
end

inputData.cfg = cfg;

disp(' ')

if cfg.writeFiles
    disp(['Writing inputData to .mat file in ',cfg.output_fd,'...'])
    cd(cfg.output_fd)
    save(cfg.output_fn,'inputData')
else
    disp('WARNING: You have selected not to save the inputData')
end

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                    End of inputData run                     ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
toc

if cfg.writeFiles; diary off; end
cd(iWasHere)