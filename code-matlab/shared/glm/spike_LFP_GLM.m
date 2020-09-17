function sd = rhythmGLMfit_shuf(cfg_in)
% GLM for spike prediction with LFP features
%
% as rhythmGLMfit, but add exclusion of duplicate cells
%
% This top-level function fits a number of GLMs to single session spike train data.
%
% The overall idea is to test whether addition of LFP features, such as
% theta power and/or gamma phase, improve cross-validated model prediction.
%
% The analysis proceeds as follows:
% - load data
% - define models (in sd.m.modelspec), MUST include a baseline model
% - prepare session-wide variables (linearized position, LFP features,
%   speed etc) on common timebase ('TVECc')
% - for each cell:
%   * prepare regressors for this ell
%   - for each cross-validation run ("pleat"; a set of folds) and fold:
%     + fit models on training data
%     + test models
%     + store error relative to baseline model
% - for each model:
%   * plot error across cells
%   * plot error across time to reward (tuning curve)
%
% run this with data session to analyze as the working folder.
%
% required path: striatal-spike-rhythms github repo

%% params
cfg_master = []; % overall params
cfg_master.dt = 0.001;
cfg_master.f.lowGamma = [40 65];
cfg_master.f.highGamma = [70 100];
cfg_master.nPleats = 1;
cfg_master.kFold = 2;
cfg_master.plotOutput = 1;
cfg_master.writeOutput = 0;
cfg_master.writeFullError = 0; % write full nSamples x nCells x nModels error for each model (NOTE: takes up to 1GB per session)
cfg_master.smooth = 501; % smoothing window (in samples) for error
cfg_master.output_prefix = 'S0_';
cfg_master.output_dir = 'C:\temp';
cfg_master.ttr_bins = [-5:0.1:5]; % time bin edges for time-to-reward tuning curves
cfg_master.linposBins = 101; % number of position bins (is autoscaled for each session)
cfg_master.nMinSpikes = 100; % minimum number of spikes needed to include cell
cfg_master.ccMethod = 'MvdM'; % cell type classification method
cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
cfg_master.iS = []; % current session number out of fd list, get this from input cfg
cfg_master.fd = []; % full list of session fd's, get this from input cfg
cfg_master.fd_extra = []; % get this from input cfg
cfg_master.nShuf = 1;

cfg_master = ProcessConfig(cfg_master,cfg_in);

%% loading
% load data
LoadExpKeys;

% spikes
sd.S = LoadSpikesTarget(cfg_master);
nSpikes = cellfun(@length, sd.S.t); keep = nSpikes >= cfg_master.nMinSpikes;
sd.S = SelectTS([], sd.S, keep);
nCells = length(sd.S.t); if nCells == 0, sd = []; return; end

%% Categorize cells and add tetrode depths
cfg_wv = []; cfg_wv.cMethod = cfg_master.ccMethod;
s_out = CategorizeStriatumWave(cfg_wv, sd.S);

s_out.unit = [s_out.other s_out.msn s_out.fsi];
s_out.ident = [zeros(1, length(s_out.other)) ones(1, length(s_out.msn)) repmat(2, 1, length(s_out.fsi))];

cfg_tt = []; cfg_tt.verbose = 1;
cfg_tt.this_rat_ID = cfg_master.fd_extra.ratID_num(cfg_master.iS);
cfg_tt.this_date = cfg_master.fd_extra.fd_date_num(cfg_master.iS);

for iC = 1:length(sd.S.t)
    sd.S.usr.cell_type(iC) = s_out.ident(find(s_out.unit == iC));
    sd.S.usr.tetrodeDepths(iC) = ExpKeys.TetrodeDepths(sd.S.usr.tt_num(iC));
    
    cfg_tt.ttno = sd.S.usr.tt_num(iC);
    [sd.S.usr.distanceTurned(iC), prev_fd] = DistanceTurned(cfg_tt, cfg_master.fd, cfg_master.fd_extra);
    cfg_tt.verbose = 0;
end

% correlate with previous session waveforms if available
if isempty(prev_fd) % no previous day available
    sd.S.usr.duplicate = zeros(size(sd.S.usr.tt_num));
else
    pushdir(prev_fd);
    S2 = LoadSpikes([]);
    nSpikes = cellfun(@length, S2.t); keep = nSpikes >= cfg_master.nMinSpikes;
    S2 = SelectTS([], S2, keep);
    
    s_out2 = CategorizeStriatumWave(cfg_wv, S2);
    s_out = CalcWVDistances([], s_out, s_out2); % add comparison with previous day's waveforms
    
    popdir;
    
    % for each cell in current session, decide if duplicate
    for iC = 1:length(sd.S.t)
        
        this_tt_no = sd.S.usr.tt_num(iC);
        prev_day_cells = find(S2.usr.tt_num == this_tt_no);
        
        if isempty(prev_day_cells) % no cells recorded fron this tt in previous session
            sd.S.usr.duplicate(iC) = 0;
        else % previous day cells found
            temp_corr = s_out.corr(iC, prev_day_cells);
            temp_peakn = s_out.peakdiffn(iC, prev_day_cells);
            
            if temp_corr > cfg_master.maxPrevCorr & abs(temp_peakn) < cfg_master.maxPeakn % wv correlation big, peak difference small
                sd.S.usr.duplicate(iC) = 1;
            else
                sd.S.usr.duplicate(iC) = 0;
            end
        end
    end
end % of previous day available checks

% LFP - vStr
if isfield(ExpKeys,'goodGamma_vStr')
    cfg = []; cfg.fc = ExpKeys.goodGamma_vStr(1);
elseif isfield(ExpKeys,'goodGamma')
     cfg = []; cfg.fc = ExpKeys.goodGamma(1);
else
    error('Don''t know what LFP to load.');
end
csc = LoadCSC(cfg);

lfp_tt = regexp(cfg.fc, 'CSC\d+', 'match');
lfp_tt = str2double(lfp_tt{1}{1}(4:end)); % need this to skip cells from same tt (could make into function)
fprintf('LFP ttno is %d\n', lfp_tt);

cfg_phi = []; % LFP features
cfg_phi.dt = median(diff(csc.tvec));
cfg_phi.ord = 100;
cfg_phi.bins = -pi:pi/18:pi;
cfg_phi.interp = 'nearest';

% restrict spikes to a timeWindow of +/-5 seconds around the rewardsite-1 
rt = getRewardTimes();
rt = rt(rt > ExpKeys.TimeOnTrack);

% Sometimes (in R117-2007-06-12, for instance) getRewardTimes() returns
% times that are spaced out less than 5 sec apart (possibly erroneus). 
% Getting rid of such reward times to maintain consistency
rt_dif = diff(rt);
rt_dif = find(rt_dif <= 5);
valid_rt = true(length(rt),1);
for i = 1:length(rt_dif)
    valid_rt(rt_dif(i)) = false;
    valid_rt(rt_dif(i)+1) = false;
end
rt = rt(valid_rt);

% For near reward_trials  
w_start = rt - 5;
w_end = rt + 5;
rt_iv = iv(w_start, w_end);
mrt_iv = MergeIV([], rt_iv);
if length(mrt_iv.tstart) ~= length(rt_iv.tstart)
    % Something funny in this session, skip for nown
    break
end

od.S = restrict(sd.S, rt_iv);

% Print number of near trials
fprintf('%d near trials found.\n', length(rt_iv.tstart));


