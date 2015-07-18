% session plotter
clear all; pack

%% load data
cfg = [];
cfg.load_questionable_cells = 1;
cfg.whichEvents = 'all';
cfg.plotOutput = 0;
cfg.removeTheta = 1;
cfg.removeRunning = 1;

load(FindFile('*metadata.mat'));
run(FindFile('*keys.m'));

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

cfg_SWR = []; cfg_SWR.fc = ExpKeys.goodSWR(1);
lfp_SWR = LoadCSC(cfg_SWR);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% linearization block %%
%%%%%%%%%%%%%%%%%%%%%%%%%
pos = LoadPos([]);

% get similar trial numbers
if strcmp(S_orig.cfg.SessionID(1:4),'R042')
    metadata = TrimTrialTimes([],metadata); % R042 only!!
end
[L_trl,R_trl] = GetMatchedTrials([],metadata,ExpKeys);
% [L_trl,R_trl] = GetMatchedTrials([],metadata);

S_left = restrict(S_orig,L_trl); posL = restrict(pos,L_trl);
S_right = restrict(S_orig,R_trl); posR = restrict(pos,R_trl);

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

clear CoordL CoordLrs CoordR CoordRrs S_left S_right posL posL_binned posR posR_binned nBins trial_ID t

%% apply speed filter to encoding data
spd = getLinSpd([],pos);
cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 10; run_iv = TSDtoIV(cfg_spd,spd);

ENC_data(1).S = restrict(ENC_data(1).S,run_iv);
ENC_data(2).S = restrict(ENC_data(2).S,run_iv);

% keep this for later display
S_run = restrict(S_orig,metadata.taskvars.trial_iv);

%% Make tuning curves for left or right trajectorie 
clear TC %in case you are re-running
cfg_tc = [];
cfg_tc.binSize = 1;

TC(1) = MakeTC(cfg_tc,ENC_data(1).S,ENC_data(1).pos);
TC(2) = MakeTC(cfg_tc,ENC_data(2).S,ENC_data(2).pos);


cfg_tc = []; cfg_tc.order = 2; cfg_tc.binsize = binSize; cfg_tc.cp = cpL;% fields (1 is peaks)

if cfg.plotOutput, TCPlot(cfg_tc,TC(1)); end

cfg_tc.cp = cpR;
if cfg.plotOutput, TCPlot(cfg_tc,TC(2)); end

% NOTE this takes a lot of Java heap space..
if cfg.plotOutput, fh = AnotherTCplotTwin([],S_run,TC(1),TC(2),pos); else fh = []; end

%% for now, only do left data...
%TC = TC_left;
%iT = 1;

%% restrict data to ordered place cells in arms only -- NOTE need a function for subsetting ts, tsd, iv etc
clear S_pc DEC_data;
for iT = 1:2

    field_order = TC(iT).field_template_idx;
    
    S_pc(iT) = S_orig;
    S_pc(iT).t = S_pc(iT).t(field_order); S_pc(iT).label = S_pc(iT).label(field_order);
    
    % restrict encoding data to improve ordering with CCFs
    ENC_data(iT).S.t = ENC_data(iT).S.t(field_order);
    ENC_data(iT).S.label = ENC_data(iT).S.label(field_order);

end

%% reorder based on CCF (during runs only) -- this doesn't work great for R050-2014-04-02 (right)
for iT = 1:2
    
    cfg_ccf = []; cfg_ccf.PlotOutput = 0; cfg_ccf.InteractiveMode = 0;
    %cfg.tch = tch;
    [~,lags1] = FieldOrderCCF(cfg_ccf,ENC_data(iT).S);
    
    %% reorder data and TCs
    ENC_data(iT).S.t = ENC_data(iT).S.t(lags1.perm_idx);
    ENC_data(iT).S.label = ENC_data(iT).S.label(lags1.perm_idx);
    
    TC(iT).field_template_idx = TC(iT).field_template_idx(lags1.perm_idx);
    TC(iT).field_loc = TC(iT).field_loc(lags1.perm_idx);
    
    S_pc(iT).t = S_pc(iT).t(lags1.perm_idx); S_pc(iT).label = S_pc(iT).label(lags1.perm_idx);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% get SWR candidates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
load(FindFile('*candidates.mat'));
if isfield(evt,'data'), evt = rmfield(evt,'data'); end

fprintf('nEvents loaded: %d\n',length(evt.tstart));

switch cfg.whichEvents
    case 'all'
        evt = restrict(evt,ExpKeys.prerecord(1),ExpKeys.postrecord(2));
    case 'prerecord'
        evt = restrict(evt,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
    case 'postrecord'
        evt = restrict(evt,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
end
fprintf('nEvents after cfg.whichEvents restrict: %d\n',length(evt.tstart));


%% remove some events
% remove events during theta
if cfg.removeTheta
    fprintf('current nEvents: %d\n',length(evt.tstart));
    
    cfg_theta = []; cfg_theta.f = [6 10]; cfg_theta.order = 4; cfg_theta.display_filter = 0; cfg_theta.type = 'fdesign';
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

%%
cfg = []; %cfg.evt = evt; 
cfg.lfp = lfp_SWR; cfg.lfpHeight = 5; cfg.lfpMax = 5;
MultiRaster(cfg,S_pc(1))

%% prepare for RUN screenshot
axis off; axis tight
set(gcf,'Color',[1 1 1])
set(gcf,'InvertHardcopy','off')

set(gca,'XLim',[8027.28 8033.35]); xl = xlim; % for R050-2014-04-02, left
hold on;
scalebar_h = plot([xl(1) xl(1)+1],[0 0],'LineWidth',2,'Color',[0.5 0.5 0.5]);

print(gcf,'-dpng','-r300','run.png');
print(gcf,'-dill','run.ai');

%% SWR screenshot
cfg = []; %cfg.evt = evt; 
cfg.lfp = lfp_SWR;  cfg.lfpHeight = 5; cfg.lfpMax = 5;
MultiRaster(cfg,S_pc(2))

axis off; axis tight
set(gcf,'Color',[1 1 1])
set(gcf,'InvertHardcopy','off')

set(gca,'XLim',[8568.8 8569.8]); xl = xlim; % for R050-2014-04-02, right
hold on;
scalebar_h = plot([xl(1) xl(1)+0.1],[0 0],'LineWidth',2,'Color',[0.5 0.5 0.5]);

print(gcf,'-dpng','-r300','swr.png');
print(gcf,'-dill','swr.ai');

%%
MultiRasterTwin(cfg,S_pc_left,S_pc_right)