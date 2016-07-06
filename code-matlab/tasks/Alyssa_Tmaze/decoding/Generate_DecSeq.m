function out = Generate_DecSeq(cfg_in)
% function out = Generate_DecSeq(cfg_in)
%
% run decoding sequence detection analysis from within session data folder
%
% assumes cfg_in.output_dir exists; image and data files are written into this
%
% CONFIGS:
%
% cfg_def = [];
% cfg_def.output_dir = 'files';
% cfg_def.output_file_prefix = []; % when writing files, prefix this
% cfg_def.dt = 0.025; % either vary this (keeping others constant), or
% cfg_def.TCsmooth = 1; % bins SD; vary this and the next
% cfg_def.QsmoothSD = 0.002; % time SD; vary this with previous
% cfg_def.minSpikes = 25; % remove cells with less than this number of spikes during run
% cfg_def.nBins = 112; % ~3cm bins for full maze (334 cm); NOTE this is number of bin edges
% cfg_def.trialStartOffset = -1; % start "run" at this time relative to center pb break (in s)
% cfg_def.tc_baseline = 0.1; % baseline firing rate, replaces zeros in TC when unsmoothed with smoothed Q
% cfg_def.load_questionable_cells = 1;
% cfg_def.includeAllCells = 1; % otherwise, only place cells
% cfg_def.trackExcludeStart = 20; % exclude this amount (in cm) from start and end of track
% cfg_def.trackExcludeEnd = 15;
% cfg_def.nSpikesHist = -0.5:105;
% cfg_def.nMinNeurons = 4;
% cfg_def.maxJump_cm = 40;
% cfg_def.minSeqLength = 10; % note, this is in bins
% cfg_def.plotOutput = 0;
% cfg_def.Qdt = 0.005;
% cfg_def.Qboxcar = 5;
% cfg_def.writeFiles = 1;
% cfg_def.removeInterneurons = 1;
%
% MvdM 2015

%% master config
cfg_def = [];
cfg_def.output_dir = 'files';
cfg_def.output_file_prefix = []; % when writing files, prefix this
cfg_def.dt = 0.025; % either vary this (keeping others constant), or
cfg_def.TCsmooth = 1; % bins SD; vary this and the next
cfg_def.QsmoothSD = 0.002; % time SD; vary this with previous
cfg_def.minSpikes = 25; % remove cells with less than this number of spikes during run
cfg_def.nBins = 112; % ~3cm bins for full maze (334 cm); NOTE this is number of bin edges
cfg_def.trialStartOffset = -1; % start "run" at this time relative to center pb break (in s)
cfg_def.tc_baseline = 0.1; % baseline firing rate, replaces zeros in TC when unsmoothed with smoothed Q
cfg_def.load_questionable_cells = 1;
cfg_def.includeAllCells = 1; % otherwise, only place cells
cfg_def.trackExcludeStart = 20; % exclude this amount (in cm) from start and end of track
cfg_def.trackExcludeEnd = 15;
cfg_def.nSpikesHist = -0.5:105;
cfg_def.nMinNeurons = 4;
cfg_def.maxJump_cm = 40;
cfg_def.minSeqLength = 10;
cfg_def.plotOutput = 0;
cfg_def.Qdt = 0.005;
cfg_def.Qboxcar = 5;
cfg_def.writeFiles = 1;
cfg_def.removeInterneurons = 1;

nMaxLaps = 20;
%cfg_def.encdecmat = 1-eye(20);
cfg_def.encdecmat = ones(1,nMaxLaps);

cfg = ProcessConfig(cfg_def,cfg_in);

%% load data
please = []; please.load_questionable_cells = cfg.load_questionable_cells;
S = LoadSpikes(please);

LoadExpKeys;
LoadMetadata;

if cfg.removeInterneurons
    cfg_temp = []; cfg_temp.fc = ExpKeys.goodSWR(1);
    csc = LoadCSC(cfg_temp);
    S = RemoveInterneuronsHC([],S,csc);
end

cfg_pos = []; cfg_pos.convFact = ExpKeys.convFact;
pos = LoadPos(cfg_pos); % pos is now in cm
        
%% set up data structs for L, R laps
clear expCond;
expCond(1).label = 'left';
expCond(2).label = 'right';

% match trials
[left,right] = GetMatchedTrials([],metadata,ExpKeys);
expCond(1).t = left;
expCond(2).t = right;

nLapsMax = max(length(expCond(1).t.tstart),length(expCond(2).t.tstart));

% tighter run boundaries: start run based on center photobeam break
evt = getEvents_Tmaze();
for iCond = 1:length(expCond)
    
    this_runStart = expCond(iCond).t.tstart;
    for iT = 1:length(this_runStart)
        next_centerpb_break_idx = nearest_idx3(this_runStart,evt.center_pb,1);
        expCond(iCond).t.tstart = evt.center_pb(next_centerpb_break_idx)'+cfg.trialStartOffset;
    end
end

expCond(1).coord = metadata.coord.coordL_cm;
expCond(2).coord = metadata.coord.coordR_cm;

expCond(1).S = S; % for making tuning curves, "encoding" model
expCond(2).S = S; % this gets restricted by rat running, on track, etc

expCond(1).decS = S; % for decoding
expCond(2).decS = S; % this only gets selections by nSpikes, etc.

%% set up output paths
this_fd = pwd;
output_fd = cat(2,pwd,'\',cfg.output_dir);
base_fn = cat(2,cfg.output_file_prefix,S.cfg.SessionID);

%% linearize paths (snap x,y position samples to nearest point on experimenter-drawn idealized track)
fprintf('Linearizing...');


chp = tsd(0,metadata.coord.chp_cm,{'x','y'}); % make choice point useable by cobebase functions

nCond = length(expCond);
for iCond = 1:nCond
    
    cfg_linpos = []; cfg_linpos.Coord = expCond(iCond).coord;
    expCond(iCond).linpos = LinearizePos(cfg_linpos,pos);
    
    % ensure that linpos is now in cm
    expCond(iCond).linpos.data = (expCond(iCond).linpos.data ./ length(cfg_linpos.Coord)).*ExpKeys.pathlength;
    
    % get cp in linpos coordinates
    expCond(iCond).cp = LinearizePos(cfg_linpos,chp);
    expCond(iCond).cp.data = (expCond(iCond).cp.data ./ length(cfg_linpos.Coord)).*ExpKeys.pathlength;
        
end
    
%% find intervals where rat is running
spd = getLinSpd([],pos); % get speed (in "camera pixels per second")

cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 5;
run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above 5 pix/s

%% exclude positions at beginning and end of track
cfg_track1 = []; cfg_track1.method = 'raw'; cfg_track1.threshold = cfg.trackExcludeStart;
cfg_track2 = []; cfg_track2.method = 'raw'; cfg_track2.dcn = '<'; cfg_track2.threshold = ExpKeys.pathlength - cfg.trackExcludeEnd;

for iCond = 1:nCond
    track_iv1 = TSDtoIV(cfg_track1,expCond(iCond).linpos);
    
    fh = @(x) restrict(x,track_iv1);
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    track_iv2 = TSDtoIV(cfg_track2,expCond(iCond).linpos);
    
    fh = @(x) restrict(x,track_iv2);
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
end

%% deal with some rat specific oddities
[~,fn,~] = fileparts(pwd);
switch fn(1:4)
    case 'R042'
        % deal with non-spike sorted chewing intervals
        tf = FindFile('*times.mat');
        if isempty(tf)
            error('No *times.mat file found for R042');
        else
            load(tf);
            chew_iv = iv(t_start*10^-6,t_end*10^-6); % times are in neuralynx timestamps, so convert to s
        end
        
    case 'R044'
        % deal with HS detachments -- remove extended times with no spikes
        cfg_Q = []; cfg_Q.dt = 1;
        Q = MakeQfromS(cfg_Q,S);
        spk_count = tsd(Q.tvec,sum(Q.data));
        
        cfg_det = []; cfg_det.threshold = 0.5; cfg_det.dcn = '>'; cfg_det.method = 'raw';
        chew_iv = TSDtoIV(cfg_det,spk_count);
        
    otherwise
        chew_iv = [];
        
end

%% restrict (linearized) position data and spike data to desired intervals
fprintf('Restricting data...');
for iCond = 1:nCond
    
    fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    if ~isempty(chew_iv)
        fh = @(x) restrict(x,chew_iv); % restrict S and linpos to non-detached times
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    end
    
    % also remove cells with insufficient spikes
    [expCond(iCond).S,cell_keep_idx] = removeEmptyCells(expCond(iCond).S);
    expCond(iCond).decS = SelectTS([],expCond(iCond).decS,cell_keep_idx);
    
    spk_count = getSpikeCount([],expCond(iCond).S);
    cell_keep_idx = spk_count >= cfg.minSpikes;
    
    expCond(iCond).S = SelectTS([],expCond(iCond).S,cell_keep_idx);
    expCond(iCond).decS = SelectTS([],expCond(iCond).decS,cell_keep_idx);
end

%% get tuning curves & decode
for iCond = 1:nCond
    
    this_TCsmooth = cfg.TCsmooth;
    
    %% construct tuning curves
    cfg_tc = [];
    cfg_tc.binEdges{1} = linspace(0,ExpKeys.pathlength,cfg.nBins); % ~3cm bins for full maze (smaller for R042...)
    
    if this_TCsmooth ~= 0
        cfg_tc.smoothingKernel = gausskernel(51,this_TCsmooth);
    end
    
    enc_laps = cfg.encdecmat(1:length(expCond(iCond).t.tstart));  % find laps to use for decoding
    cfg_select = []; cfg_select.verbose = 0;
    lap_iv = SelectIV(cfg_select,expCond(iCond).t,find(enc_laps));
    
    enc_S = restrict(expCond(iCond).S,lap_iv);
    enc_linpos = restrict(expCond(iCond).linpos,lap_iv);
    
    % could be empty -- no cells for this lap's encoding model
    if all(cellfun(@isempty,enc_S.t))
        continue;
    end
    
    expCond(iCond).tc = TuningCurves(cfg_tc,enc_S,enc_linpos);
    
    % keep track of cp
    [~,expCond(iCond).cp_bin] = histc(expCond(iCond).cp.data,cfg_tc.binEdges{1});
    
    if cfg.plotOutput
        figure;
        imagesc(expCond(iCond).tc.tc);
    end
    
    % get Q-matrix
    this_Qsd = cfg.QsmoothSD;
    
    
    %% Q-mat
    cfg_Q = [];
    cfg_Q.dt = cfg.Qdt;
    cfg_Q.boxcar_size = cfg.Qboxcar;
    
    if this_Qsd == 0
        cfg_Q.smooth = [];
    else
        cfg_Q.smooth = 'gauss';
        cfg_Q.gausswin_sd = this_Qsd;
    end
    
    expCond(iCond).Q = MakeQfromS(cfg_Q,expCond(iCond).decS);
    
    % raw Q-mat, to get nNeurons later
    cfg_Qraw = []; cfg_Qraw.smooth = []; cfg_Qraw.dt = cfg_Q.dt; cfg_Qraw.boxcar_size = cfg_Q.boxcar_size;
    expCond(iCond).Qraw = MakeQfromS(cfg_Qraw,expCond(iCond).decS);
    expCond(iCond).nNeurons = sum(expCond(iCond).Qraw.data >= 1);
    
    %% decode
    cfg_decode = [];
    cfg_decode.nMinSpikes = cfg.dt;
    cfg_decode.excludeMethod = 'frate';
    expCond(iCond).P = DecodeZ(cfg_decode,expCond(iCond).Q,expCond(iCond).tc.tc);
    
    %% quantify decoding accuracy on RUN
    this_trueZ = tsd(expCond(iCond).linpos.tvec,expCond(iCond).tc.pos_idx);
    cfg_err = []; cfg_err.mode = 'max';
    
    keep_idx = unique(nearest_idx3(this_trueZ.tvec,expCond(iCond).P.tvec));
    this_Pscore = expCond(iCond).P;
    this_Pscore.tvec = this_Pscore.tvec(keep_idx);
    this_Pscore.data = this_Pscore.data(:,keep_idx);
    expCond(iCond).Perr = DecodeErrorZ(cfg_err,this_Pscore,this_trueZ);
    
    %% exclude bins with insufficient neurons
    expCond(iCond).P2 = expCond(iCond).P;
    expCond(iCond).P2.data(:,expCond(iCond).nNeurons < cfg.nMinNeurons) = NaN;
    
    %% sequence detection
    cfg_seq = [];
    binSize_cm = ExpKeys.pathlength/cfg.nBins; % cm in one bin
    cfg_seq.maxJump = round(cfg.maxJump_cm/binSize_cm); % number of bins corresponding to 40cm
    cfg_seq.minLength = cfg.minSeqLength;
    [expCond(iCond).seq_iv,expCond(iCond).decode_map] = DecSeqDetectZ(cfg_seq,expCond(iCond).P2);
    
    %% plot
    if cfg.plotOutput
        figure;
        
        cfg_plot = [];
        cfg_plot.bgcol = '.k';
        cfg_plot.fgcol = '.r';
        a1 = subplot(311);
        PlotTSDfromIV(cfg_plot,expCond(iCond).seq_iv,expCond(iCond).decode_map);
        
        a2 = subplot(312);
        imagesc(expCond(iCond).Qraw.tvec,1:size(expCond(iCond).Qraw.data,1),expCond(iCond).Qraw.data);
        
        a3 = subplot(313);
        plot(expCond(iCond).Qraw.tvec,expCond(iCond).nNeurons,'.');
        
        linkaxes([a1 a2 a3],'x');
    end
    
    %% remove source matrices (too memory-intensive)
    % NOTE: could keep Q and decoded bits for detected sequences
    % (helpful for shuffle checks later)?
    expCond(iCond).Q = []; expCond(iCond).Qraw = []; expCond(iCond).P = []; expCond(iCond).P2 = [];
    
end % of iCond loop

out.expCond = expCond;
out.fn = fn;
out.cfg = cfg;

%% write data -- could be helper function
if cfg.writeFiles
    cd(output_fd);
    out_fn = cat(2,base_fn,'-DecSeq_data.mat');
    save(out_fn,'out');
    cd(this_fd)
end