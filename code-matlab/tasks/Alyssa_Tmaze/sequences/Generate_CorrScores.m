function out = Generate_CorrScores(cfg_in)
% function out = Generate_CorrScores(cfg_in)
%
% automated workflow for correlation-based sequence scoring, run from within data folder
%
% assumes cfg_in.output_dir exists; image files are written into this
%
% CONFIGS:
%
% cfg_def.load_questionable_cells = 1;
% cfg_def.output_dir = 'files'; % subdir to write files to
% cfg_def.output_file_prefix = []; % prefix use when writing files
% cfg_def.whichEvents = 'all'; % 'prerecord', 'postrecord', 'all'
% cfg_def.dt = 0.025; % time step of starts to consider within candidate iv
%  window
% cfg_def.twin = 0.025:0.025:0.2; % window sizes (from start to end) to
%  consider within candidate iv window
% cfg_def.nActiveCells = 4; % need MORE than this number of cells active
%  (>, not >=)
% cfg_def.writeFiles = 1;
% cfg_def.plotOutput = 0;

cfg_def = [];
cfg_def.load_questionable_cells = 1;
cfg_def.output_dir = 'files'; % subdir to write files to
cfg_def.output_file_prefix = 'CS_'; % prefix use when writing files
cfg_def.whichEvents = 'all'; % 'prerecord', 'postrecord', 'all'
cfg_def.dt = 0.025; % time step of starts to consider within candidate iv window
cfg_def.twin = 0.025:0.025:0.2; % window sizes (from start to end) to consider within candidate iv window
cfg_def.nActiveCells = 4; % need MORE than this number of cells active (>, not >=)
cfg_def.writeFiles = 1;
cfg_def.plotOutput = 0;
cfg_def.removeTheta = 1;
cfg_def.removeRunning = 1;
cfg_def.matchFields = 0;

cfg = ProcessConfig2(cfg_def,cfg_in);

%
LoadMetadata
LoadExpKeys

% Load spikes
please = [];
please.useClustersFile = 0;
please.load_questionable_cells = cfg.load_questionable_cells;
please.getTTnumbers = 1;
S_orig = LoadSpikes(please);

% load theta CSC for exclusion (note could do this once and store the
% times...)
cfg_theta = []; cfg_theta.fc = ExpKeys.goodTheta(1);
lfp_theta = LoadCSC(cfg_theta);

% set up filenames and paths
this_fd = pwd;
output_fd = cat(2,pwd,'\',cfg.output_dir);
base_fn = cat(2,cfg.output_file_prefix,S_orig.cfg.SessionID);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% linearization block %%
%%%%%%%%%%%%%%%%%%%%%%%%%
pos = LoadPos([]);

% get similar trial numbers
if strcmp(S_orig.cfg.SessionID(1:4),'R042')
    metadata = TrimTrialTimes([],metadata); % R042 only!!
end
[L_trl,R_trl] = GetMatchedTrials([],metadata,ExpKeys);
% [L_trl,R_trl] = GetMatchedTrials_old([],metadata);

S_left = restrict2(S_orig,L_trl); posL = restrict2(pos,L_trl);
S_right = restrict2(S_orig,R_trl); posR = restrict2(pos,R_trl);

% linearize runs and bin positions
CoordL = metadata.coord.coordL; CoordR = metadata.coord.coordR;

% check if coord actually overlaps with real position data
if cfg.plotOutput
    figure(1); subplot(221);
    plot(getd(posL,'x'),getd(posL,'y'),'.k'); hold on;
    plot(CoordL(1,:),CoordL(2,:),'og'); title('LeftCoord');
    
    subplot(222);
    plot(getd(posR,'x'),getd(posR,'y'),'.k'); hold on;
    plot(CoordR(1,:),CoordR(2,:),'og'); title('RightCoord');
end

% standardize Coord to have specific bin size
run_dist = ExpKeys.pathlength; % distance travelled on a single run of the track in cm (T-maze)
binSize = 3; % in cm (2 for D&G, 1 for Davidson et al)
nBins = round(run_dist/binSize);
CoordLrs(1,:) = interp1(1:size(CoordL,2),CoordL(1,:),linspace(1,size(CoordL,2),nBins),'linear');
CoordLrs(2,:) = interp1(1:size(CoordL,2),CoordL(2,:),linspace(1,size(CoordL,2),nBins),'linear');
CoordRrs(1,:) = interp1(1:size(CoordR,2),CoordR(1,:),linspace(1,size(CoordR,2),nBins),'linear');
CoordRrs(2,:) = interp1(1:size(CoordR,2),CoordR(2,:),linspace(1,size(CoordR,2),nBins),'linear');

cfg_c = []; cfg_c.Coord = CoordLrs;
posL_binned = LinearizePos(cfg_c,posL);

% cp
cpL = tsd(0,metadata.coord.chp,{'x','y'});
cpL = LinearizePos(cfg_c,cpL); cpL = cpL.data(1);

cfg_c = []; cfg_c.Coord = CoordRrs;
posR_binned = LinearizePos(cfg_c,posR);

% cp
cpR = tsd(0,metadata.coord.chp,{'x','y'});
cpR = LinearizePos(cfg_c,cpR); cpR = cpR.data(1);

% store left/right linearized data in a single struct (cleaner workspace!)
ENC_data(1).trial_type = 'left';
ENC_data(1).Coord = CoordLrs;
ENC_data(1).pos = posL_binned;
ENC_data(1).S = S_left;

ENC_data(2).trial_type = 'right';
ENC_data(2).Coord = CoordRrs;
ENC_data(2).pos = posR_binned;
ENC_data(2).S = S_right;

% check binned position (should be no/few gaps)
if cfg.plotOutput
    figure(1);
    subplot(223);
    plot(ENC_data(1).pos.tvec,getd(ENC_data(1).pos,'z'),'.')
    title('Left-linearized pos');
    subplot(224);
    plot(ENC_data(2).pos.tvec,getd(ENC_data(2).pos,'z'),'.')
    title('Right-linearized pos');
end

% write fig
if cfg.writeFiles, cd(output_fd);
    out_fn = cat(2,base_fn,'-CoordCheck.png');
    print(gcf,'-r300','-dpng',out_fn);
    close(gcf);
    cd(this_fd);
end

clear CoordL CoordLrs CoordR CoordRrs S_left S_right posL posL_binned posR posR_binned nBins trial_ID t

%% apply speed filter to encoding data
spd = getLinSpd([],pos);
cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 10; run_iv = TSDtoIV(cfg_spd,spd);

ENC_data(1).S = restrict2(ENC_data(1).S,run_iv);
ENC_data(2).S = restrict2(ENC_data(2).S,run_iv);

% keep this for later display
S_run = restrict2(S_orig,metadata.taskvars.trial_iv);

%% Make tuning curves for left or right trajectorie 
clear TC %in case you are re-running
cfg_tc = [];
cfg_tc.binSize = 1;

TC(1) = MakeTC(cfg_tc,ENC_data(1).S,ENC_data(1).pos);
TC(2) = MakeTC(cfg_tc,ENC_data(2).S,ENC_data(2).pos);
TC(1).tc = TC(1).tc';
TC(2).tc = TC(2).tc';

if cfg.matchFields % NOTE, could consider this after going to fields beyond CP only...
   [TC(1),TC(2)] = GetMatchedFields([],TC(1),TC(2));
end

cfg_tc = []; cfg_tc.order = 2; cfg_tc.binsize = binSize; cfg_tc.cp = cpL;% fields (1 is peaks)

if cfg.plotOutput, TCPlot(cfg_tc,TC(1)); end

% write fig -- could be helper function
if cfg.writeFiles, cd(output_fd);
    out_fn = cat(2,base_fn,'-TC_left.png');
    print(gcf,'-r300','-dpng',out_fn);
    close(gcf);
    cd(this_fd);
end

cfg_tc.cp = cpR;
if cfg.plotOutput, TCPlot(cfg_tc,TC(2)); end
% write fig -- could be helper function
if cfg.writeFiles, cd(output_fd);
    out_fn = cat(2,base_fn,'-TC_right.png');
    print(gcf,'-r300','-dpng',out_fn);
    close(gcf);
    cd(this_fd);
end

% NOTE this takes a lot of Java heap space..
if cfg.plotOutput, fh = AnotherTCplotTwin([],S_run,TC(1),TC(2),pos); else fh = []; end
% write figs
if cfg.writeFiles, cd(output_fd);
    for iFH = 1:length(fh)
        out_fn = cat(2,base_fn,'-fields_',num2str(iFH),'.png');
        set(fh(iFH), 'InvertHardCopy', 'off');
        print(fh(iFH),'-r300','-dpng',out_fn);
        close(fh(iFH));
    end
    cd(this_fd);
end

%% for now, only do left data...
%TC = TC_left;
%iT = 1;

%% restrict data to ordered place cells in arms only -- NOTE need a function for subsetting ts, tsd, iv etc
clear S_pc DEC_data;
for iT = 1:2

    field_order = TC(iT).field_template_idx(TC(iT).field_loc > cpL);
    
    S_pc(iT) = S_orig;
    S_pc(iT).t = S_pc(iT).t(field_order); S_pc(iT).label = S_pc(iT).label(field_order);
    
    % restrict encoding data to improve ordering with CCFs
    ENC_data(iT).S.t = ENC_data(iT).S.t(field_order);
    ENC_data(iT).S.label = ENC_data(iT).S.label(field_order);
    
    % decode with place cell data only
    DEC_data(iT).S = S_orig;
    DEC_data(iT).S.t = DEC_data(iT).S.t(field_order);
    DEC_data(iT).S.label = DEC_data(iT).S.label(field_order);

end

%% reorder based on CCF (during runs only) -- this doesn't work great for R050-2014-04-02 (right)
for iT = 1:2
    
    cfg_ccf = []; cfg_ccf.PlotOutput = 0; cfg_ccf.InteractiveMode = 0;
    %cfg.tch = tch;
    [~,lags1] = FieldOrderCCF(cfg_ccf,ENC_data(iT).S);
    
    %% reorder data and TCs
    ENC_data(iT).S.t = ENC_data(iT).S.t(lags1.perm_idx);
    ENC_data(iT).S.label = ENC_data(iT).S.label(lags1.perm_idx);
    
    DEC_data(iT).S.t = DEC_data(iT).S.t(lags1.perm_idx);
    DEC_data(iT).S.label = DEC_data(iT).S.label(lags1.perm_idx);
    
    TC(iT).template_idx = TC(iT).template_idx(lags1.perm_idx);
    
    TC(iT).peak_idx = TC(iT).peak_idx(lags1.perm_idx);
    TC(iT).peak_loc = TC(iT).peak_loc(lags1.perm_idx);
    TC(iT).tc = TC(iT).tc(:,TC(iT).template_idx);
    TC(iT).spk_hist = TC(iT).spk_hist(TC(iT).template_idx,:);
    
    S_pc(iT).t = S_pc(iT).t(lags1.perm_idx); S_pc(iT).label = S_pc(iT).label(lags1.perm_idx);
    
    out.ENC_data = ENC_data;
end
%% How to plot reordered tuning cirves? should be option in TCPlot...

%% Extract SWR events (or load)
%load(FindFile('*candidates.mat')); % loads an evt structure with candidates
LoadCandidates
evt = rmfield(evt,'data');

fprintf('nEvents loaded: %d\n',length(evt.tstart));

switch cfg.whichEvents
    case 'all'
        %evt = restrict(evt,ExpKeys.prerecord(1),ExpKeys.postrecord(2));
    case 'prerecord'
        evt = restrict(evt,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
    case 'postrecord'
        evt = restrict(evt,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
    case 'task'
        evt = restrict(evt,ExpKeys.task(1),ExpKeys.task(2));
end
fprintf('nEvents after cfg.whichEvents restrict: %d\n',length(evt.tstart));


%% remove some events
% remove events during theta
if cfg.removeTheta
    fprintf('current nEvents: %d\n',length(evt.tstart));
    
    cfg_theta = []; cfg_theta.f = [6 10]; cfg_theta.order = 4; cfg_theta.display_filter = 0; %cfg_theta.type = 'fdesign';
    lfp_theta = FilterLFP(cfg_theta,lfp_theta);
    
    tpow = LFPpower([],lfp_theta);
    
    cfg_theta = []; cfg_theta.method = 'zscore'; cfg_theta.threshold = 2; cfg_theta.dcn =  '<';
    low_theta_iv = TSDtoIV(cfg_theta,tpow);
    
    evt = restrict(evt,low_theta_iv);
    fprintf('nEvents after theta filtering: %d\n',length(evt.tstart));
end

% remove events during running
if cfg.removeRunning
    fprintf('current nEvents: %d\n',length(evt.tstart));
    
    cfg_run = []; cfg_run.method = 'raw'; cfg_run.threshold = 10; cfg_run.dcn = '<';
    run_iv = TSDtoIV(cfg_run,spd);
    
    evt = restrict(evt,run_iv);
    fprintf('nEvents after speed filtering: %d\n',length(evt.tstart));
end

%% add some padding
cfg_pad = []; cfg_pad.d = [-0.1 0.1];
evt = ResizeIV(cfg_pad,evt);

%% remove events with less than some number of cells
cfg_n = []; cfg_n.label = 'NActiveCells_left';
evt = AddNActiveCellsIV(cfg_n,evt,S_pc(1)); % note this has binning at 5ms by default
cfg_n = []; cfg_n.label = 'NActiveCells_right';
evt = AddNActiveCellsIV(cfg_n,evt,S_pc(2)); % note this has binning at 5ms by default

keep_idx = find(evt.usr.('NActiveCells_left') > cfg.nActiveCells | evt.usr.('NActiveCells_right') > cfg.nActiveCells);
evt = SelectIV([],evt,keep_idx);
fprintf('nEvents after nActiveCells filtering: %d\n',length(evt.tstart));

%% plot SWR events
% iT = 1;
% cfg = [];
% % cfg.lfp = lfp;
% cfg.lfp = lfp_theta;
% cfg.evt = evt_lt;
% cfg.evt = manual_iv;
% cfg.windowSize = 0.5; %in seconds
% MultiRaster(cfg,S_pc(1)); %NOTE: MultiRaster has the navigate function inside!!!

% % or..
% cfg = [];
% cfg.evt = evt;
% cfg.lfp = lfp_theta;
% PolyRaster(cfg,S_pc(1),S_pc(2));

%% Run sequence analysis

%make sure that there are any significant sequences
if isempty(evt.tstart)
    out.score1 = nan;
    out.score3 = nan;
    out.sign_evt = nan;
else
    for iT = 1:2
        cfg_cs_shuf3 = [];
        cfg_cs_shuf3.dt = cfg.dt;
        cfg_cs_shuf3.twin = cfg.twin;
        cfg_cs_shuf3.shuffleType = 3;
        %     cfg_cs_shuf3.nShuffles = 10;
        out.score3(iT) = CorrScoreWin3(cfg_cs_shuf3,evt,S_pc(iT));
        
        cfg_cs_shuf1 = cfg_cs_shuf3;
        cfg_cs_shuf1.shuffleType = 1;
        %     cfg_cs_shuf1.nShuffles = 10;
        out.score1(iT) = CorrScoreWin3(cfg_cs_shuf1,evt,S_pc(iT));
        
        keep_idx = find(out.score1(iT).WIN_rho_perc > 0.95 & out.score3(iT).WIN_rho_perc > 0.95);
        %     keep_idx = find(score1(iT).WIN_rho_perc < 0.05 & score3(iT).WIN_rho_perc < 0.05);
        
        sign_evt = iv;
        sign_evt.tstart = vertcat(out.score1(iT).WIN_iv(keep_idx).tstart);
        sign_evt.tend = vertcat(out.score1(iT).WIN_iv(keep_idx).tend);
        sign_evt.usr = out.score1(iT).WIN_iv(1).usr;
        sign_evt.cfg = out.score1(iT).WIN_iv(1).cfg;
        
        % sign_evt = evt_lt;
        % sign_evt.tstart = sign_evt.tstart(keep_idx);
        % sign_evt.tend = sign_evt.tend(keep_idx);
        %
        % sign_evt = rmfield(sign_evt,'usr');
        % sign_evt.usr(1).label = 'rho (obs)';
        % sign_evt.usr(1).data = score1(iT).WIN_rho_obs(keep_idx);
        % sign_evt.usr(2).label = 'idshuf perc';
        % sign_evt.usr(2).data = 1-score1(iT).WIN_rho_perc(keep_idx);
        % sign_evt.usr(3).label = 'tshuf perc';
        % sign_evt.usr(3).data = 1-score3(iT).WIN_rho_perc(keep_idx);
        
        % cfg_SWR = []; cfg_SWR.fc = ExpKeys.goodSWR(1);
        % lfp_SWR = LoadCSC(cfg_SWR);
        
        % cfg = [];
        % cfg.lfp = lfp_SWR;
        % cfg.evt = sign_evt;
        % cfg.windowSize = 0.5; %in seconds
        % MultiRaster(cfg,S_pc(iT)); %NOTE: MultiRaster has the navigate function inside!!!
    end
    
end

out.cfg = cfg; 
out.nLcells = length(S_pc(1).t); 
out.nRcells = length(S_pc(2).t);
out.cand_evt = evt;

%% write data -- could be helper function
if cfg.writeFiles
    cd(output_fd);
    out_fn = cat(2,base_fn,'-CorrScores.mat');
    save(out_fn,'out');
    
    if ~isempty(evt.tstart)
        Plot_SingleCorrScore;
        out_fn = cat(2,base_fn,'-corrSummary.png');
        set(gcf, 'InvertHardCopy', 'off');
        print(gcf,'-r300','-dpng',out_fn);
        close all;
    end
    cd(this_fd)
end
