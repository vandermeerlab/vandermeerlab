%% sandbox to generate cross correlations
%
%
% The analysis proceeds as follows:
% - load data
% - remove spikes arouhnd +/-5 sec window around reward times and retain
% neurons with at least 100 spikes after the removal
% - separate type 1 (MSNs) and type 2 (FSI) neurons
% - generate cross correlation MSN pairs (cx1), FSI pairs (cx2), and
%    MSN-FSI pairs (cx3)
%
% output variables are stored in the od variable, which can be saved and
% later used to plot the figures
% 
% TODO: Use restict() and invertIV() to remove spikes around the reward
% times

%% setup
clear;
cd('/Users/manishm/Work/vanDerMeerLab/ADRLabData/');
please = [];
please.rats = {'R117'};%,'R119','R131','R132'}; % vStr-only rats
[cfg_in.fd,cfg_in.fd_extra] = getDataPath(please);
cfg_in.write_output = 1;
cfg_in.output_dir = '/Users/manishm/Work/vanDerMeerLab/Common/temp';
cfg_in.cx_binsize = 0.01; 

%%
% Top level loop which calls the main function for all the sessions
for iS = 1:2%length(cfg_in.fd) % for each session...
    
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
    cfg_master.cx_binsize = 0.01;

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

    %restrict spikes to times that animal was on the track
    sd.S = restrict(sd.S, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);

    % reward deliveries
    evt = LoadEvents([]);
    reward_t = evt.t{1}; % should generalize this with known labels etc.. and remove double labels

    % Remove spikes 5 sec from around the reward times
    for i = 1:length(sd.S.t)
        sd.S.t{i} = RemoveSpikes(sd.S.t{i}, reward_t, 5);
    end

    sd.S.cell_type = sd.S.usr.cell_type;

    exc_types = 0;
    %Keep cells greater with greater than 100 spikes and of the allowed types
    sd.S = KeepCells(sd.S,cfg_master.nMinSpikes,exc_types);

    %Crosscorrelation for MSNs
    c1 = sd.S.t(sd.S.cell_type == 1);
    c2 = sd.S.t(sd.S.cell_type == 2);
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

    % 10 millsecond bins for ccf
    cfg_cx.binsize = cfg_master.cx_binsize;

    k = 1;
    for i = 1:(n1-1)
        for j = (i+1):n1
            [od.cx1{k},~] = ccf(cfg_cx,c1{i},c1{j});
            k = k+1;
        end
    end

    k = 1;
    for i = 1:(n2-1)
        for j = (i+1):n2
            [od.cx2{k},~] = ccf(cfg_cx,c2{i},c2{j});
            k = k+1;
        end
    end

    k = 1;
    for i = 1:n1
        for j = 1:n2
            [od.cx3{k},~] = ccf(cfg_cx,c1{i},c2{j});
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

%%
% function to remove spikes from ost that are in a window of length w_length
% around timpepoints in r_time
% very inefficient implementation
function rst = RemoveSpikes(ost, r_time, w_length) 
    keep = true(length(ost), 1);
    for iS = 1:length(ost)
            cur_spk = ost(iS);
            for iT = 1:length(r_time)
                low = r_time(iT) - w_length;
                high = r_time(iT) + w_length;
                if cur_spk >= low && cur_spk <= high
                    keep(iS)= false;
                    break
                end
            end        
    end
    rst = ost(keep);
end
