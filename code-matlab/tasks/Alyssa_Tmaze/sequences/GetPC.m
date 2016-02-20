function out = GetPC(cfg_in)
% function out = GetPC(cfg_in)
%
% get S_pc(1) (left) and S_pc(2) (right)
%
% CONFIGS:
%
% cfg_def.load_questionable_cells = 1;

cfg_def = [];
cfg_def.load_questionable_cells = 1;
cfg_def.matchFields = 0;

cfg = ProcessConfig2(cfg_def,cfg_in);

%
load(FindFile('*metadata.mat'));
run(FindFile('*keys.m'));

% Load spikes
please = [];
please.useClustersFile = 0;
please.load_questionable_cells = cfg.load_questionable_cells;
please.getTTnumbers = 1;
S_orig = LoadSpikes(please);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% linearization block %%
%%%%%%%%%%%%%%%%%%%%%%%%%
pos = LoadPos([]);

% get similar trial numbers
if strcmp(S_orig.cfg.SessionID(1:4),'R042')
    metadata = TrimTrialTimes([],metadata); % R042 only!!
end
%[L_trl,R_trl] = GetMatchedTrials([],metadata,ExpKeys); % AC version - errors?
[L_trl,R_trl] = GetMatchedTrials_old([],metadata); % MvdM old version that includes bad trials...

S_left = restrict(S_orig,L_trl); posL = restrict(pos,L_trl);
S_right = restrict(S_orig,R_trl); posR = restrict(pos,R_trl);

% linearize runs and bin positions
CoordL = metadata.coord.coordL; CoordR = metadata.coord.coordR;

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

TC(1) = MakeTC(cfg_tc,ENC_data(1).S,ENC_data(1).pos); TC(1).tc = TC(1).tc';
TC(2) = MakeTC(cfg_tc,ENC_data(2).S,ENC_data(2).pos); TC(2).tc = TC(2).tc';


if cfg.matchFields % NOTE, could consider this after going to fields beyond CP only...
   [TC(1),TC(2)] = GetMatchedFields([],TC(1),TC(2));
end

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
    
end

out.S_pc = S_pc;