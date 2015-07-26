function [data_out] = SimSWRstats_all(cfg_in,data_in)
%% get SWR event data
cfg_def = [];
cfg_def.requireExpKeys = 1;
cfg_def.ExpKeysFields = {'goodSWR','goodTheta'};
cfg_def.requireMetadata = 1;
cfg_def.MetadataFields = {'SWRfreqs'};
cfg_def.requireVT = 1;
cfg_def.plotOutput = 0;
cfg_def.rats = {'R042','R044','R050'};
cfg_def.restrict = [];
cfg_def.NAU = 0;
originalFolder = pwd;

cfg = ProcessConfig2(cfg_def,cfg_in);


%% Generate simulated data and grab statistics
for iR = 1:length(cfg.rats)
    rn = cfg.rats{iR};
    cfg_temp.rats = cfg.rats(iR);
    fd = sort(getTmazeDataPath(cfg_temp)); % get all session directories  
    
    simdata.TC.(rn) = [];
    simdata.ENC_data.(rn) = [];
    simdata.R(2).(rn) = [];
    simdata.SWR_size(2).(rn) = [];
    simdata.IFFR(2).(rn) = [];
    simdata.IFFR_field(2).(rn) = [];
    simdata.IFFR_fieldmax(2).(rn) = [];
    simdata.ISI_SWR(2).(rn) = [];
    simdata.ISI_IF(2).(rn) = [];
    simdata.path_len(2).(rn) = [];
    simdata.nCells = [];
    simdata.nSWR = [];
    
    %% run for each session
    for iFD = 1:length(fd)
        dirname = fd{iFD}
        cd(dirname);

        LoadExpKeys
        LoadMetadata % for freqs

        %load candidate events
        load(FindFile('*candidates.mat'));
        if strcmp(cfg.restrict,'pre');
            evt = restrict2(evt,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
        elseif strcmp(cfg.restrict,'post');
            evt = restrict2(evt,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
        elseif strcmp(cfg.restrict,'task');
            evt = restrict2(evt,ExpKeys.task(1),ExpKeys.task(2));
        end
        
        %load empirical data
        ENC_data = data_in.ENC_data.(rn){1,iFD};
        cfg_S = [];
        cfg_S.useClustersFile = 0;
        cfg_S.load_questionable_cells = 0;
        cfg_S.min_cluster_quality = 3;
        cfg_S.getTTnumbers = 1;
        S_orig = LoadSpikes(cfg_S);
        
        TC_temp = data_in.TC.(rn){1,iFD};
        
        % get similar trial numbers
        if strcmp(S_orig.cfg.SessionID(1:4),'R042')
            metadata = TrimTrialTimes([],metadata); % R042 only!!
        end
        [L_trl,R_trl] = GetMatchedTrials([],metadata,ExpKeys);
        
        simdata.TC.(rn){iFD} = TC_temp;
        simdata.ENC_data.(rn){iFD} = ENC_data;
        
        %% Generate stats for each arm
        run_trl = [L_trl R_trl];
        for iT=1:2
            %% Get SWR specific stats
            % restrict data to place cells
            field_order = TC_temp(iT).template_idx;
            if isempty(field_order)
                sprintf('No place cells!')
                continue 
            elseif length(field_order) < 2
                sprintf('Not enough place cells!')
                continue
            end
            
            S_pc(iT) = S_orig;
            S_pc(iT).t = S_pc(iT).t(field_order); 

            iv_in = iv(evt.tstart,evt.tend);
            if isempty(iv_in.tstart)
                sprintf('no candidates!')
                continue
            end
            
            if cfg.NAU
                % limit candidates by NAU of place cells
                minCells = 4;
                activeCellsIV = AddNActiveCellsIV([],iv_in,S_pc(iT));
                evt.nActiveCells = activeCellsIV.usr.data;
                exclude = arrayfun(@(x) x < minCells,evt.nActiveCells);
                evt.tstart(exclude) = [];
                evt.tend(exclude) = [];
                evt.nActiveCells(exclude) = [];
            end

            nCells = length(field_order)
            nSWR = length(evt.tstart)
            
            if nSWR == 0
                sprintf('no significant events!')
                continue
            end
            
            % make position vector from runs
            pos_in = tsd();
            pos_in.tvec = ENC_data(iT).pos.tvec;
            pos_in.data = getd(ENC_data(iT).pos,'z');
            pos_in = restrict2(pos_in,run_trl(iT).tstart(1),run_trl(iT).tend(1));
            
            % make place fields
            for i=1:length(field_order)
                curr_field = TC_temp(iT).tc(:,field_order(i));
                pfs(i).xc = TC_temp(iT).peak_loc{i};          
                pfs(i).maxFiringRate = curr_field(pfs(i).xc);
                pfs(i).rates = curr_field;
                assert(length(pfs(i).xc) == length(pfs(i).maxFiringRate),'!');
            end
            
            R_temp = zeros(length(field_order),nSWR);
            SWR_size_temp = zeros(nSWR,1);
            ISI_temp = [];
            path_length_temp = [];
            
            for iTrial = 1:nSWR %for each trial...
                % get empirical "path length" for each SWR
                pathlen = data_in.path_len(iT).(rn){1,iFD}(iTrial);
                if pathlen<1 || isnan(pathlen)
                    continue
                end
                
                % generate virtual path based on a segment of a real run compressed by SWR length
                path_idx = pos_in.data >= max(pos_in.data)-pathlen; % distance from reward
                tlen = data_in.SWR_size(iT).(rn){1,iFD}(iTrial);
                pos_virtual = pos_in;
                pos_virtual.data = pos_in.data(path_idx);
                pos_virtual.tvec = linspace(0,tlen,length(pos_virtual.data));

                assert(~isempty(pos_virtual.data),'empty virtual run')
                
                % simulate place cell firing for event
                cfg_model.nTrials = 1;
                cfg_model.lockingParameter = 0; % no theta modulation
                cfg_model.convertfile = 'spikes';
                cfg_model.verbose = 0;
                cfg_model.rate_dt = 0.0001;
                S_temp = pfmodel2(cfg_model,pfs,pos_virtual);

                % gather data across cells
                spk_t = [];
                for iC=1:length(S_temp.t) %for each cell...
                    R_temp(iC,iTrial) = length(S_temp.t{iC});
                    spk_t = vertcat(spk_t,S_temp.t{iC});
                    ISI_temp = vertcat(ISI_temp,diff(S_temp.t{iC}));
                end %iterate cells
                                
                % get length of event
                if isempty(spk_t) || length(spk_t) < 2
                    SWR_size_temp(iTrial) = nan;
                else
                    SWR_size_temp(iTrial) = max(spk_t) - min(spk_t);
                end
                
                % get virtual path length of each event (% max dist between first and last peak)
                path_inds = find(R_temp(:,iTrial));
                if isempty(path_inds);
                    path_length_temp = vertcat(path_length_temp,nan); 
                    continue
                end
                peak_loc_temp = [];
                for iP = 1:length(path_inds) % for each cell find peak
                    peak_loc_temp = horzcat(peak_loc_temp,pfs(iP).xc);
                end
                path_start = min(peak_loc_temp); 
                path_end = max(peak_loc_temp);
                path_length_temp = vertcat(path_length_temp,(path_end-path_start)); 
                
            end %iterate trials

            simdata.R(iT).(rn){iFD} = R_temp;
            simdata.SWR_size(iT).(rn){iFD} = SWR_size_temp;
            simdata.ISI_SWR(iT).(rn){iFD} = ISI_temp;
            simdata.path_len(iT).(rn){iFD} = path_length_temp;
            simdata.nCells = vertcat(simdata.nCells,nCells);
            simdata.nSWR = vertcat(simdata.nSWR,nSWR);
            
            clear R_temp S_temp pfs
        end %iterate left vs right
        
    end %iterate session
    
end %iterate rat

data_out = simdata;
cd(originalFolder)

end