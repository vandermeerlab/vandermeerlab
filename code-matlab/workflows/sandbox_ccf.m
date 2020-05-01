%% sandbox to generate cross correlations
%
%
% The analysis proceeds as follows:
% - load data
% - restrict spikes to around +/-5 sec window around reward times and retain
%   neurons with at least 100 spikes remaining
% - separate type 1 (MSNs) and type 2 (FSI) neurons
% - generate cross correlation MSN pairs (cx1), FSI pairs (cx2), and
%    MSN-FSI pairs (cx3)
%
% output variables are stored in the od variable, which can be saved and
% later used to plot the figures

%% setup
clear;
cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/');
please = [];
please.rats = {'R117','R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.output_dir = '/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp';
cfg_in.cx_binsize = 0.01; 

%%
% Top level loop which calls the main function for all the sessions
for iS = 1:length(cfg_in.fd) % for each session...
    
    cfg_in.iS = iS;
    pushdir(cfg_in.fd{iS});
    generateCCF(cfg_in); % do the business
    popdir;
    
end % of sessions

%%
% Main function to generate cross-correlations
function od = generateCCF(cfg_in)

    LoadExpKeys;
    
    % params
    cfg_master = []; % overall params
    cfg_master.dt = 0.001;
    cfg_master.ccMethod = 'MvdM'; % cell type classification method
    cfg_master.maxPrevCorr = 0.99; % if wv correlation with previous day is bigger than this, cell is possible duplicate
    cfg_master.maxPeakn = 0.2; % if peak wv difference (normalized) with previous day is smaller than this, cell is possible duplicate
    cfg_master.nMinSpikes = 100;
    cfg_master.iS = 1; % current session number out of fd list, get this from input cfg
    cfg_master.fd = []; % full list of session fd's, get this from input cfg
    cfg_master.fd_extra = []; % get this from input cfg
    cfg_master.write_output = 0;
    cfg_master.output_prefix = 'ccf_';
    cfg_master.output_dir = 'C:\temp';
    cfg_master.exc_types = 0; %cell types to be excluded
    cfg_master.cx_msn = 0.01;  % bin size for cross-correlation between pairs of Medium Spiny Neurons (MSNs) in seconds
    cfg_master.cx_fsi = 0.005;  % bin size for cross-correlations between pairs of Fast Spiking Interneurons (FSIs) in seconds
    cfg_master.cx_mix = 0.005;  % bin size for cross-correlations between pairs of 1 FSI abd 1 MSN in seconds
    cfg_master.max_t = 0.5; % half window length for the cross correlations in seconds in seconds

    cfg_master = ProcessConfig(cfg_master,cfg_in);
    
    % spikes
    sd.S = LoadSpikesTarget(cfg_master);

    % Categorize cells and add tetrode depths
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

    % restrict spikes to a timeWindow of +/-5 seconds around the reward
    % times
    rt = getRewardTimes();
    w_start = rt - 5;
    w_end = rt + 5;
    rt_iv = iv(w_start, w_end);
    rt_iv = MergeIV([],rt_iv);
    sd.S = restrict(sd.S, rt_iv);

    sd.S.cell_type = sd.S.usr.cell_type;
    sd.S.tt_id = sd.S.usr.tt_num;

    exc_types = 0;
    % Keep cells greater with greater than nMinSpike spikes and of the allowed types
    sd.S = KeepCells(sd.S,cfg_master.nMinSpikes,exc_types);
    
    % Setting up parameters for cross_correlations
    c1 = sd.S.t(sd.S.cell_type == 1);
    c2 = sd.S.t(sd.S.cell_type == 2);
    od.tt1 = sd.S.tt_id(sd.S.cell_type == 1);
    od.tt2 = sd.S.tt_id(sd.S.cell_type == 2);
    n1 = length(c1);
    t1 = (n1*(n1-1))/2;
    n2 = length(c2);
    t2 = (n2*(n2-1))/2;
    n3 = n1*n2;
    od.cx1 = cell(t1,1);
    od.cx2 = cell(t2,1);
    od.cx3 = cell(n3,1);
    od.l1 = sd.S.label(sd.S.cell_type == 1);
    od.l2 = sd.S.label(sd.S.cell_type == 2);
    od.tvec1 = (-cfg_master.max_t:cfg_master.cx_msn:cfg_master.max_t);
    od.tvec2 = (-cfg_master.max_t:cfg_master.cx_fsi:cfg_master.max_t);
    od.tvec3 = (-cfg_master.max_t:cfg_master.cx_mix:cfg_master.max_t);
    
    cfg_cx.max_t = cfg_master.max_t;
    cfg_cx.smooth = 1;
    % Cross correlations for MSNs
    cfg_cx.binsize =  cfg_master.cx_msn;
    cfg_cx.gauss_w = 7*cfg_cx.binsize;
    cfg_cx.gauss_sd = cfg_cx.binsize;
    k = 1;
    for i = 1:(n1-1)
        for j = (i+1):n1
            if od.tt1(i) == od.tt1(j) %same tetrode
                od.cx1{k} = nan(length(od.tvec1),1);
            else
                [od.cx1{k},~] = ccf2(cfg_cx,c1{i},c1{j});
            end
            k = k+1;
        end
    end
    
    % Cross correlations for FSIs
    cfg_cx.binsize =  cfg_master.cx_fsi;
    cfg_cx.gauss_w = 15*cfg_cx.binsize;
    cfg_cx.gauss_sd = cfg_cx.binsize;
    k = 1;
    for i = 1:(n2-1)
        for j = (i+1):n2
            if od.tt2(i) == od.tt2(j) %same tetrode
                od.cx2{k} = nan(length(od.tvec2),1);
            else
            [od.cx2{k},~] = ccf2(cfg_cx,c2{i},c2{j});
            end
            k = k+1;
        end
    end

    % Cross correlations for MSN-FSI pairs
    cfg_cx.binsize =  cfg_master.cx_mix;
    cfg_cx.gauss_w = 15*cfg_cx.binsize;
    cfg_cx.gauss_sd = cfg_cx.binsize;
    k = 1;
    for i = 1:n1
        for j = 1:n2
            if od.tt1(i) == od.tt2(j) %same tetrode
                od.cx3{k} = nan(length(od.tvec3),1);
            else
            [od.cx3{k},~] = ccf2(cfg_cx,c1{i},c2{j});
            end
            k = k+1;
        end
    end
    if cfg_master.write_output
         [~, fp, ~] = fileparts(pwd);
         pushdir(cfg_master.output_dir);
         fn_out = cat(2,cfg_master.output_prefix, fp, '_od.mat');
         save(fn_out,'od'); % should add option to save in specified output dir
         popdir;
    end
end

% Other functions
%%
function S = KeepCells(S,minSpikes,cTypes)
    thr = minSpikes;

    nCells = length(S.t);
    for i = nCells:-1:1
       l(i) = length(S.t{i}); 
    end
    keep = l >= thr;
    for iC = 1:length(keep)
        for iT = 1:length(cTypes)
            if S.cell_type(iC) == cTypes(iT)
                keep(iC) = 0;
                break
            end
        end
    end
    S.label = S.label(keep);
    S.t = S.t(keep);
    S.cell_type = S.cell_type(keep);
    S.tt_id = S.tt_id(keep);
end


%%
function S = LoadSpikesTarget(cfg_in)
    if ~isfield(cfg_in, 'Target') % no target specified, load them all
        S = LoadSpikes([]);
        return;
    end
    LoadExpKeys;
    % target specified, need to do some work
    % first see if this session has more than one target
    nTargets = length(ExpKeys.Target);
    if ~iscell(ExpKeys.Target) || (iscell(ExpKeys.Target) && length(ExpKeys.Target) == 1) % one target
        target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
        if isempty(target_idx)
            S = ts;
        else
            S = LoadSpikes([]);
        end
    else % multiple targets, assume TetrodeTargets exists
        please = []; please.getTTnumbers = 1;
        S = LoadSpikes(please);
        target_idx = strmatch(cfg_in.Target, ExpKeys.Target);
        tt_num_keep = find(ExpKeys.TetrodeTargets == target_idx);
        keep = ismember(S.usr.tt_num, tt_num_keep);
        S = SelectTS([], S, keep);
    end
end
