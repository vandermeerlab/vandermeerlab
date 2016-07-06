%%
clear all;

cfg = [];
cfg.prefix = 'R2_'; % which files to load
cfg.whichEvents = 'postrecord'; % {'prerecord','taskrest','taskrun','postrecord'}; % which events to process?
cfg.whichSeq = 'either'; % {'all','fwd','bwd','either'}; % which sequences to process?
cfg.sessions = {'food','water'};
cfg.arms = {'left','right'};
cfg.writeDiary = 1; % keep a text file record of command window history
cfg.output_fd = 'C:\projects\AlyssaTmaze\resultsFiles'; % where to place output files
cfg.saveData = 1;
cfg.rats = {'R042','R044','R050','R064'};
cfg.cpbin = 10; % if non-empty, restrict to sequences with start or end beyond this (relative to CP)
cfg.minActiveCells = 5;
cfg.minlen = 0.05; % otherwise, minimum length in s
cfg.output_prefix = cat(2,cfg.prefix,'DecSeq_',cfg.whichEvents,'_',cfg.whichSeq,'_');
cfg.SWRoverlap = 1; % if 1, only keep events detected as SWR candidates
%cfg.output_prefix = cat(2,cfg.prefix,'DecSeq_',cfg.whichEvents,'_',cfg.whichSeq,'_cp70_');


%%
fd = getTmazeDataPath(cfg);
iWasHere = pwd;

if cfg.writeDiary % save command window text
    warning off
    cd(cfg.output_fd)
    diary([cfg.output_prefix,'stats','.txt'])
    cd(cfg.output_fd)
    disp(' ')
    disp(date)
    disp(' ')
end

disp(['Script name: ',mfilename])
disp(' ')
disp('You have selected: ')
disp(cfg.whichEvents)
disp(cfg.whichSeq)
disp(cfg.rats)
disp(' ')

%% GET COMBINED DATA
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%                                                               %%%')
disp('%%%               Combined sequences data (all rats):             %%%')
disp('%%%                                                               %%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% init vars to place data into
ALL_sig_seq = [];

ALL_sig_seq.count = []; % number of sig events
ALL_sig_seq.countN = []; % proportions of sig events
ALL_sig_seq.arm = []; % left (1), right (2)
ALL_sig_seq.type = []; % restriction type (food/water)
ALL_sig_seq.sess = []; % session ID in case we want to restrict later
ALL_sig_seq.decErr = []; % track decoding error during RUN
ALL_sig_seq.firstChoice = []; % left (1), right (2)

ALL_sig_seq.rsq = [];
ALL_sig_seq.pval = [];
ALL_sig_seq.beta = [];
ALL_sig_seq.len = [];

swapfun = @(x) -x+3; % utility function to access variable from other condition


%% collect
nSessions = 0;
for iFD = 1:length(fd)
    nSessions = nSessions +1;
    close all;
    
    %%% load data for this session %%%
    cd(fd{iFD});
    LoadExpKeys;
    LoadCandidates;
    LoadMetadata;
    
    %cd('files');
    cd([fd{iFD},'\files'])
    [~,sessionID,~] = fileparts(fd{iFD}); % 'RXXX-201X-XX-XX'
    this_file = FindFiles([cfg.prefix,sessionID,'-DecSeq_data.mat']);
    
    if isempty(this_file)
        fprintf('Session %s: no DecSeq file found, skipping...\n',fd{iFD});
        continue;
    end
    
    load(this_file{1}); cd .. % loads 'out' variable saved by ALL_Generate_DecSeq.m
    
    % track some useful things about this session
    this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
    switch metadata.taskvars.sequence{1}
        case 'L'
            firstChoice = 1;
        case 'R'
            firstChoice = 2;
    end
    seq_tosave = []; % for exporting as "candidates" for later use
    
    % for each arm, count number of sequences that meet criteria
    for iCond = 1:length(cfg.arms)
        
        this_seqR = [];
        
        % only include events that overlap with SWRs
        if cfg.SWRoverlap
            this_seq = out.expCond(iCond).seq_iv;
            this_seqR = IntersectIV([],this_seq,evt);
            else
                this_seqR = out.expCond(iCond).seq_iv;
        end
        
        % only select events with some minimum length
        if ~isempty(cfg.minlen)
            evt_len = this_seqR.tend-this_seqR.tstart;
            this_seqR = SelectIV([],this_seqR,evt_len >= cfg.minlen);
        end
        
        seq_tosave = UnionIV([],seq_tosave,this_seqR); % for exporting as "candidates", keep both L and R
        
        % exclude events that are sequences for both
        other_seq = out.expCond(swapfun(iCond)).seq_iv;
        other_seqR = IntersectIV([],other_seq,evt);
        
        this_seqR = DifferenceIV([],this_seqR,other_seqR);
        
        % restrict events to epoch of interest
        switch cfg.whichEvents
            case 'prerecord'
                this_seqR = restrict(this_seqR,0,ExpKeys.prerecord(2));
            case 'postrecord'
                this_seqR = restrict(this_seqR,ExpKeys.postrecord(1),Inf);
            case 'taskrest'
                this_seqR = restrict(this_seqR,metadata.taskvars.rest_iv);
            case 'taskrun'
                this_seqR = restrict(this_seqR,metadata.taskvars.trial_iv);
        end
        
        % add fwd/bwd information
        this_map = out.expCond(iCond).decode_map;
        for iS = 1:length(this_seqR.tstart)
            
            this_seq = restrict(this_map,this_seqR.tstart(iS),this_seqR.tend(iS));
            
            [betas,~,~,~,stat] = regress(this_seq.data',[ones(1,length(this_seq.data)); 1:length(this_seq.data)]');
            
            this_seqR.usr.beta(iS) = betas(2); % regression coefficient
            this_seqR.usr.rsq(iS) = stat(1); % variance explained
            this_seqR.usr.pval(iS) = stat(3); % p-value
            
            if ~isempty(cfg.cpbin) % does this sequence include a piece beyond the choice point?
                this_seqR.usr.loc(iS) = double(this_seq.data(1) > out.expCond(iCond).cp_bin + cfg.cpbin | this_seq.data(end) > out.expCond(iCond).cp_bin + cfg.cpbin); % start or end should be beyond cp
            end
            
        end
        
        if ~isempty(cfg.cpbin)  % restrict to having at least some piece on L/R arms
            cfg_s = []; cfg_s.operation = '>'; cfg_s.threshold = 0.5;
            this_seqR = SelectIV(cfg_s,this_seqR,'loc');
        end
        
        if ~isempty(this_seqR.tstart)
            switch cfg.whichSeq % select fwd/bwd if requested
                
                case 'all'
                    keep_idx = 1:length(this_seqR.usr.pval);
                case 'fwd'
                    keep_idx = (this_seqR.usr.pval < 0.05 & this_seqR.usr.beta > 0);
                case 'bwd'
                    keep_idx = (this_seqR.usr.pval < 0.05 & this_seqR.usr.beta < 0);
                case 'either'
                    keep_idx = (this_seqR.usr.pval < 0.05);
                    
            end
            this_seqR = SelectIV([],this_seqR,keep_idx);
        end
        
        % only select events with some minimum number of cells active
        % note this is in addition to decoding-level selection which
        % requires some number of neurons active in a bin
        % (nMinNeurons in Generate_DecSeq)
        
        if ~isempty(this_seqR.tstart)
            this_seqR = AddNActiveCellsIV([],this_seqR,out.expCond(iCond).decS);
            
            keep_idx = this_seqR.usr.nActiveCells >= cfg.minActiveCells;
            this_seqR = SelectIV([],this_seqR,keep_idx);
        end
        
        % count events
        nEvents = length(this_seqR.tstart);
        
        
        % track decoding error
        this_decErr = nanmean(out.expCond(iCond).Perr.data);
        
        ALL_sig_seq.decErr = cat(1,ALL_sig_seq.decErr,this_decErr);
        ALL_sig_seq.count = cat(1,ALL_sig_seq.count,nEvents);
        ALL_sig_seq.countN = cat(1,ALL_sig_seq.countN,nEvents);
        ALL_sig_seq.arm = cat(1,ALL_sig_seq.arm,iCond);
        ALL_sig_seq.type = cat(1,ALL_sig_seq.type,this_session_type);
        ALL_sig_seq.sess = cat(1,ALL_sig_seq.sess,iFD);
        ALL_sig_seq.firstChoice = cat(1,ALL_sig_seq.firstChoice,firstChoice);
        
        % track event durations, R^2, beta, pval -- NOTE will be different
        % variable lengths than above
        if ~isempty(this_seqR.tstart)
            ALL_sig_seq.rsq = cat(1,ALL_sig_seq.rsq,this_seqR.usr.rsq');
            ALL_sig_seq.pval = cat(1,ALL_sig_seq.pval,this_seqR.usr.pval');
            ALL_sig_seq.beta = cat(1,ALL_sig_seq.beta,this_seqR.usr.beta');
            ALL_sig_seq.len = cat(1,ALL_sig_seq.len,this_seqR.tend-this_seqR.tstart);
        end
        
    end
    
    % convert raw count to proportions (left/right)
    ALL_sig_seq.countN(end-1:end) = ALL_sig_seq.countN(end-1:end)./sum(ALL_sig_seq.countN(end-1:end));
    
    % process candidates & save
    %evt = MergeIV([],seq_tosave);
    %S = LoadSpikes([]);
    %out_fn = cat(2,S.cfg.SessionID,'-DecSeqCand.mat');
    %save(out_fn,'evt');
    
end

%%
food_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 1);
food_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 1);
food_leftN = ALL_sig_seq.countN(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 1);
food_rightN = ALL_sig_seq.countN(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 1);

water_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 2);
water_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 2);
water_leftN = ALL_sig_seq.countN(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 2);
water_rightN = ALL_sig_seq.countN(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 2);

d = [nansum(food_left) nansum(food_right) nansum(water_left) nansum(water_right)];

data.all.food_left = nansum(food_left); data.all.food_leftN = nanmean(food_leftN);
data.all.food_right = nansum(food_right); data.all.food_rightN = nanmean(food_rightN);
data.all.water_left = nansum(water_left); data.all.water_leftN = nanmean(water_leftN);
data.all.water_right = nansum(water_right); data.all.water_rightN = nanmean(water_rightN);

data.all.ALL_sig_seq = ALL_sig_seq; % keep for later analysis

disp(' ')
disp('Displaying number of significant sequences: ')
disp(['Left, food restr, ',num2str(data.all.food_left)])
disp(['Right, food restr ',num2str(data.all.food_right)])
disp(['Left, water restr ',num2str(data.all.water_left)])
disp(['Right, water restr ',num2str(data.all.water_right)])
disp(' ')

% p_foodseq = tmaze_coin_pval(max([d(1) d(2)]),d(1)+d(2));
% p_waterseq = tmaze_coin_pval(max([d(3) d(4)]),d(3)+d(4));
% disp(' ');
% disp('***************************')
% disp('** binomial test results **')
% disp('***************************')
% disp(['Food day pval: ',num2str(p_foodseq)])
% disp(['Water day pval: ',num2str(p_waterseq)])

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
    ALL_sig_seq = [];
    
    ALL_sig_seq.count = []; % number of sig events
    ALL_sig_seq.countN = []; % proportions of sig events
    ALL_sig_seq.arm = []; % left (1), right (2)
    ALL_sig_seq.type = []; % restriction type (food/water)
    ALL_sig_seq.sess = []; % session ID in case we want to restrict later
    ALL_sig_seq.decErr = []; % track decoding error during RUN
    ALL_sig_seq.firstChoice = []; % left (1), right (2)
    
    ALL_sig_seq.rsq = [];
    ALL_sig_seq.pval = [];
    ALL_sig_seq.beta = [];
    ALL_sig_seq.len = [];
    
    swapfun = @(x) -x+3; % utility function to access variable from other condition
    
    
    %% collect
    nSessions = 0;
    for iFD = 1:length(fd)
        nSessions = nSessions +1;
        close all;
        
        %%% load data for this session %%%
        cd(fd{iFD});
        LoadExpKeys;
        LoadCandidates;
        LoadMetadata;
        
        %cd('files');
        cd([fd{iFD},'\files'])
        [~,sessionID,~] = fileparts(fd{iFD}); % 'RXXX-201X-XX-XX'
        this_file = FindFiles([cfg.prefix,sessionID,'-DecSeq_data.mat']);
        
        if isempty(this_file)
            fprintf('Session %s: no DecSeq_data file found, skipping...\n',fd{iFD});
            continue;
        end
        
        load(this_file{1}); cd .. % loads 'out' variable saved by ALL_Generate_DecSeq.m
        
        % track some useful things about this session
        this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
        switch metadata.taskvars.sequence{1}
            case 'L'
                firstChoice = 1;
            case 'R'
                firstChoice = 2;
        end
        seq_tosave = []; % for exporting as "candidates" for later use
        
        % for each arm, count number of sequences that meet criteria
        for iCond = 1:length(cfg.arms)
            
            this_seqR = [];
            
            % only include events that overlap with SWRs
            if cfg.SWRoverlap
                this_seq = out.expCond(iCond).seq_iv;
                this_seqR = IntersectIV([],this_seq,evt);
            else
                this_seqR = out.expCond(iCond).seq_iv;
            end
            
            % only select events with some minimum length
            if ~isempty(cfg.minlen)
                evt_len = this_seqR.tend-this_seqR.tstart;
                this_seqR = SelectIV([],this_seqR,evt_len >= cfg.minlen);
            end
            
            seq_tosave = UnionIV([],seq_tosave,this_seqR); % for exporting as "candidates", keep both L and R
            
            % exclude events that are sequences for both
            other_seq = out.expCond(swapfun(iCond)).seq_iv;
            other_seqR = IntersectIV([],other_seq,evt);
            
            this_seqR = DifferenceIV([],this_seqR,other_seqR);
            
            % restrict events to epoch of interest
            switch cfg.whichEvents
                case 'prerecord'
                    this_seqR = restrict(this_seqR,0,ExpKeys.prerecord(2));
                case 'postrecord'
                    this_seqR = restrict(this_seqR,ExpKeys.postrecord(1),Inf);
                case 'taskrest'
                    this_seqR = restrict(this_seqR,metadata.taskvars.rest_iv);
                case 'taskrun'
                    this_seqR = restrict(this_seqR,metadata.taskvars.trial_iv);
            end
            
            % add fwd/bwd information
            this_map = out.expCond(iCond).decode_map;
            for iS = 1:length(this_seqR.tstart)
                
                this_seq = restrict(this_map,this_seqR.tstart(iS),this_seqR.tend(iS));
                
                [betas,~,~,~,stat] = regress(this_seq.data',[ones(1,length(this_seq.data)); 1:length(this_seq.data)]');
                
                this_seqR.usr.beta(iS) = betas(2); % regression coefficient
                this_seqR.usr.rsq(iS) = stat(1); % variance explained
                this_seqR.usr.pval(iS) = stat(3); % p-value
                
                if ~isempty(cfg.cpbin) % does this sequence include a piece beyond the choice point?
                    this_seqR.usr.loc(iS) = double(this_seq.data(1) > out.expCond(iCond).cp_bin + cfg.cpbin | this_seq.data(end) > out.expCond(iCond).cp_bin + cfg.cpbin); % start or end should be beyond cp
                end
                
            end
            
            if ~isempty(cfg.cpbin)  % restrict to having at least some piece on L/R arms
                cfg_s = []; cfg_s.operation = '>'; cfg_s.threshold = 0.5;
                this_seqR = SelectIV(cfg_s,this_seqR,'loc');
            end
            
            if ~isempty(this_seqR.tstart)
                switch cfg.whichSeq % select fwd/bwd if requested
                    
                    case 'all'
                        keep_idx = 1:length(this_seqR.usr.pval);
                    case 'fwd'
                        keep_idx = (this_seqR.usr.pval < 0.05 & this_seqR.usr.beta > 0);
                    case 'bwd'
                        keep_idx = (this_seqR.usr.pval < 0.05 & this_seqR.usr.beta < 0);
                    case 'either'
                        keep_idx = (this_seqR.usr.pval < 0.05); 
                end
                this_seqR = SelectIV([],this_seqR,keep_idx);
            end
            
            % only select events with some minimum number of cells active
            % note this is in addition to decoding-level selection which
            % requires some number of neurons active in a bin
            % (nMinNeurons in Generate_DecSeq)
            
            if ~isempty(this_seqR.tstart)
                this_seqR = AddNActiveCellsIV([],this_seqR,out.expCond(iCond).decS);
                
                keep_idx = this_seqR.usr.nActiveCells >= cfg.minActiveCells;
                this_seqR = SelectIV([],this_seqR,keep_idx);
            end
            
            % count events
            nEvents = length(this_seqR.tstart);
            
            
            % track decoding error
            this_decErr = nanmean(out.expCond(iCond).Perr.data);
            
            ALL_sig_seq.decErr = cat(1,ALL_sig_seq.decErr,this_decErr);
            ALL_sig_seq.count = cat(1,ALL_sig_seq.count,nEvents);
            ALL_sig_seq.countN = cat(1,ALL_sig_seq.countN,nEvents);
            ALL_sig_seq.arm = cat(1,ALL_sig_seq.arm,iCond);
            ALL_sig_seq.type = cat(1,ALL_sig_seq.type,this_session_type);
            ALL_sig_seq.sess = cat(1,ALL_sig_seq.sess,iFD);
            ALL_sig_seq.firstChoice = cat(1,ALL_sig_seq.firstChoice,firstChoice);
            
            % track event durations, R^2, beta, pval -- NOTE will be different
            % variable lengths than above
            if ~isempty(this_seqR.tstart)
                ALL_sig_seq.rsq = cat(1,ALL_sig_seq.rsq,this_seqR.usr.rsq');
                ALL_sig_seq.pval = cat(1,ALL_sig_seq.pval,this_seqR.usr.pval');
                ALL_sig_seq.beta = cat(1,ALL_sig_seq.beta,this_seqR.usr.beta');
                ALL_sig_seq.len = cat(1,ALL_sig_seq.len,this_seqR.tend-this_seqR.tstart);
            end
            
        end
        
        % convert raw count to proportions (left/right)
        ALL_sig_seq.countN(end-1:end) = ALL_sig_seq.countN(end-1:end)./sum(ALL_sig_seq.countN(end-1:end));
        
        % process candidates & save
        %evt = MergeIV([],seq_tosave);
        %S = LoadSpikes([]);
        %out_fn = cat(2,S.cfg.SessionID,'-DecSeqCand.mat');
        %save(out_fn,'evt');
        
    end
    
    %%
    food_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 1);
    food_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 1);
    food_leftN = ALL_sig_seq.countN(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 1);
    food_rightN = ALL_sig_seq.countN(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 1);
    
    water_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 2);
    water_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 2);
    water_leftN = ALL_sig_seq.countN(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 2);
    water_rightN = ALL_sig_seq.countN(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 2);
    
    d = [nansum(food_left) nansum(food_right) nansum(water_left) nansum(water_right)];
    
    data.(rats{iRat}).food_left = nansum(food_left); data.(rats{iRat}).food_leftN = nanmean(food_leftN);
    data.(rats{iRat}).food_right = nansum(food_right); data.(rats{iRat}).food_rightN = nanmean(food_rightN);
    data.(rats{iRat}).water_left = nansum(water_left); data.(rats{iRat}).water_leftN = nanmean(water_leftN);
    data.(rats{iRat}).water_right = nansum(water_right); data.(rats{iRat}).water_rightN = nanmean(water_rightN);
    
    data.(rats{iRat}).ALL_sig_seq = ALL_sig_seq; % keep for later analysis
    
    disp(' ')
    disp('Displaying number of significant sequences: ')
    disp(['Left, food restr, ',num2str(data.(rats{iRat}).food_left)])
    disp(['Right, food restr ',num2str(data.(rats{iRat}).food_right)])
    disp(['Left, water restr ',num2str(data.(rats{iRat}).water_left)])
    disp(['Right, water restr ',num2str(data.(rats{iRat}).water_right)])
    disp(' ')
    
    % p_foodseq = tmaze_coin_pval(max([d(1) d(2)]),d(1)+d(2));
    % p_waterseq = tmaze_coin_pval(max([d(3) d(4)]),d(3)+d(4));
    % disp(' ');
    % disp('***************************')
    % disp('** binomial test results **')
    % disp('***************************')
    % disp(['Food day pval: ',num2str(p_foodseq)])
    % disp(['Water day pval: ',num2str(p_waterseq)])
    
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
end % of rats

%%
data.date = datestr(now);
data.script = mfilename;

if cfg.writeDiary; diary off; end

if cfg.saveData; cd(cfg.output_fd); save([cfg.output_prefix,'out.mat'],'data'); end

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

cd(iWasHere)
warning on