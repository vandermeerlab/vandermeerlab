function CAND_iv_out = getCandSeq(cfg_in,S,pos)
%% GETCANDSEQ Get Candidate Sequences
%    CAND_iv_out = getCandSeq(cfg_in,S,pos) returns an IV struct containing the time
%    intervals of candidate replay event sequences.
%
%   INPUT:
%       cfg_in - input cfg parameters
%       S - input spiketrain
%       pos - input position vector binned in cm
%
%   OUTPUT:
%       CAND_iv_out - candidate replay events as an iv struct 
%
%   CFG OPTIONS:
%       cfg.MethodSelection = [1 1 0 0];
%       cfg.spd_thr = 10; % speed threshold in cm/s%
%       cfg.SWR = swr_iv_in; %input SWR iv
%       cfg.tvec = lfp.tvec_in %if using MUA or NAU, need lfp time vector
%       cfg.mua_thr = 3;
%       cfg.mua_max_thr = 5;
%       cfg.nau_thr = 3; % >3 cells for threshold (at least 4)
%       cfg.overlap_thr = 0.01;
%       cfg.intersectFlag = 1;
%
%   Method Selection: a vector with each item being zero or one. Selection can be
%   completely restrictive (hard) or can be inclusive (soft).
%       MethodSelection(1) : Animal Speed (hard)
%       MethodSelection(2) : Multi-unit activity (soft)
%       MethodSelection(3) : Number of active units (hard)
%       MethodSelection(4) : Dragoi & Tonegawa methods (hard)
%
% youkitan 2014-12-28

%% Parse cfg inputs
cfg_def.MethodSelection = [1 1 1 0];
cfg_def.spd_thr = 10; % speed threshold in cm/s
cfg_def.mua_thr = 3;
cfg_def.mua_max_thr = 5;
cfg_def.nau_thr = 4;
cfg_def.overlap_thr = 0.01;
cfg_def.intersectFlag = 1;

cfg = ProcessConfig2(cfg_def,cfg_in);

%% Parse methods
methods = cfg.MethodSelection;

useSpd = methods(1);
useMUA = methods(2);
useNAU = methods(3);
useDragoi = methods(4);

% Create a cell array to contain all IVs
useIV = {};

%% Restrict spikes by animals running speed
if useSpd
    S_orig = S;
    vvec = getLinSpd([],pos);
    cfg_temp = [];
    cfg_temp.threshold = cfg.spd_thr;
    cfg_temp.method = 'raw';
    cfg_temp.dcn = '<';
    spd_iv = TSDtoIV(cfg_temp,vvec);
    S = restrict(S,spd_iv);
end

%% Restrict events by SWR epochs
if isfield(cfg,'SWR') && ~useMUA
    useIV = [useIV, cfg.SWR];
end

%% Restrict events by MUA
if useMUA
   if ~isfield(cfg,'tvec')
       error('Missing time vector. Add lfp time vector to cfg file! (e.g., cfg.tvec = lfp.tvec)')
   end
    
   % Setup temp cfg
    cfg_temp = [];
    cfg_temp.threshold = cfg.mua_thr; %in std, determines size of ivs
    cfg_temp.method = 'raw';
    cfg_temp.merge_thr = 0.03;
    cfg_temp.tvec = cfg.tvec;
    
    % obtain multi-unit activity and z-score
    mua = getMUA(cfg_temp,S); 
    mua_z = zscore_tsd(mua);

    % threshold mua
    muaIV = TSDtoIV(cfg_temp,mua_z);

    % to each event, add a field with the max z-scored power
    cfg_temp = [];
    cfg_temp.method = 'max'; % 'min', 'mean'
    cfg_temp.label = 'maxz'; % what to call this in iv, i.e. usr.label
    muaIV = AddTSDtoIV(cfg_temp,muaIV,mua_z);

    % select only those events of >5 z-scored power
    cfg_temp.threshold = cfg.mua_max_thr;
    muaIV = SelectIV(cfg_temp,muaIV);
    
    % select only those time intervals with high MUA that coincide with SWR events
    if isfield(cfg,'SWR')
        spwr_iv = cfg.SWR;
        numIV = length(muaIV.tstart);
        keep_idx = [];
        
        for iIV = 1:numIV
            curr_tstart = muaIV.tstart(iIV) - cfg.overlap_thr;
            curr_tend = muaIV.tend(iIV) + cfg.overlap_thr;

            check1 = any(find(curr_tstart < spwr_iv.tstart & spwr_iv.tstart < curr_tend));
            check2 = any(find(curr_tstart < spwr_iv.tend & spwr_iv.tend < curr_tend));

            if check1 || check2
                keep_idx = cat(1,keep_idx,iIV);
            end

        end

        keep_tstart = muaIV.tstart(keep_idx);
        keep_tend = muaIV.tend(keep_idx);

        muaIV = iv(keep_tstart,keep_tend);
    end
    
    useIV = [useIV, muaIV];
end

%% Restrict events by number of active units
if useNAU
    % Concatenate all spikes
    spk = [];
    for iC = 1:length(S.t)
       spk = cat(1,spk,S.t{iC}); 
    end

    % create internal tvec with 50ms bins
    bin = 0.05; %50 ms
    tveci = min(spk)-bin:bin:max(spk)+bin;
    
   % Setup temp cfg
    cfg_temp.threshold = cfg.nau_thr-1; % number of cells, determines size of thresholded regions
    cfg_temp.method = 'raw';
    cfg_temp.merge_thr = 0.05;
    cfg_temp.dt = bin;
    cfg_temp.tvec = tveci'; %must be Nx1 vector
    
    % obtain number of active units
    nau = getNAU(cfg_temp,S);
end

%% Restrict using Dragoi & Tonegawa methods
if useDragoi
    % Concatenate all spikes
    spk = [];
    for iC = 1:length(S.t)
       spk = cat(1,spk,S.t{iC}); 
    end

    % create tvec with 50ms bins
    bin = 0.05; %100 ms
    edges = min(spk):bin:max(spk);
    centers = edges(1:end-1)+bin/2;

    % count spikes per bin
    spk_binned = zeros(length(S.t),length(centers));
    for iC = 1:length(S.t)
        spk_hist = histc(S.t{iC},edges)';
        spk_hist(end-1) = spk_hist(end-1)+spk_hist(end); % combine last bin and edge
        spk_hist = spk_hist(1:end-1);
        spk_binned(iC,:) = spk_hist;
    end

    % get bins with >4 different cells firing and flanked by silent bins 
    spk_numcells = sum(spk_binned ~= 0); %number of cells firing in bin
    keep_idx = [];
    for iB = 2:length(spk_numcells)-1;
        if spk_numcells(iB) >= 4
            if spk_numcells(iB-1) ==0 && spk_numcells(iB+1) == 0
                keep_idx = [keep_idx iB];
            end 
        end
    end

    % convert to data timescale
    if ~isfield(cfg,'tvec')
       error('Missing time vector. Add lfp time vector to cfg file! (e.g., cfg.tvec = lfp.tvec)')
    end
    temp_yvals = zeros(1,length(centers));
    temp_yvals(keep_idx) = 1;
    spk_keep = interp1(centers,temp_yvals,cfg.tvec);

    % convert to IV
    spk_events_tsd = tsd;
    spk_events_tsd.tvec = cfg.tvec;
    spk_events_tsd.data = spk_keep';

    cfg_temp = [];
    cfg_temp.method = 'raw';
    cfg_temp.threshold = 0.999;
    cfg_temp.dcn =  '>';
    spk_events_iv = TSDtoIV(cfg_temp,spk_events_tsd);
    useIV = [useIV,spk_events_iv];
end


%% Create a single IV from useIV
assert(~isempty(useIV),'no interval data created!');

if length(useIV) == 1
    CAND_iv_out = useIV{1};
else
    cfg_temp = [];
    cfg_temp.intersectFlag = cfg.intersectFlag;
    CAND_iv = MergeIV(cfg_temp,useIV{:});
    CAND_iv_out = CAND_iv;
end


if useNAU
    % to each candidate event, add a field with the max number of cells firing within a window
    cfg_temp = [];
    cfg_temp.method = 'max'; % 'min', 'mean'
    cfg_temp.label = 'nCells'; % what to call this in iv, i.e. usr.label
    nauCAND_iv = AddTSDtoIV(cfg_temp,CAND_iv_out,nau);
    
    % select only those events where >nau_thr different cells fire during window
    cfg_temp.threshold = cfg.nau_thr;
    CAND_iv_out = SelectIV(cfg_temp,nauCAND_iv);
end
