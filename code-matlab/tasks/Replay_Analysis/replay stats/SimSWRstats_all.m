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
originalFolder = pwd;

cfg = ProcessConfig2(cfg_def,cfg_in);


%% Generate simulated data and grab statistics
simdata.tlen = [];
simdata.mean = [];
simdata.var = [];
simdata.prop = [];
simdata.avgspk = [];
simdata.ISI = [];


profile on
for iR = 1:length(cfg.rats)
    rn = cfg.rats{iR};
    cfg_temp.rats = cfg.rats(iR);
    fd = sort(getTmazeDataPath(cfg_temp)); % get all session directories  
    
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
        
        %% Generate stats for each arm
        run_trl = [L_trl R_trl];
        for iT=1:2
            %% Get SWR specific stats
            % restrict data to place cells
            field_order = TC_temp(iT).template_idx;
            if isempty(field_order)
                sprintf('No place cells!')
                continue 
            end
            S_pc(iT) = S_orig;
            S_pc(iT).t = S_pc(iT).t(field_order); 

            iv_in = iv(evt.tstart,evt.tend);
            if isempty(iv_in.tstart)
                sprintf('no candidates!')
                continue
            end
            
            % limit candidates by NAU of place cells
            minCells = 4;
            activeCellsIV = AddNActiveCellsIV([],iv_in,S_pc(iT));
            evt.nActiveCells = activeCellsIV.usr.data;
            exclude = arrayfun(@(x) x < minCells,evt.nActiveCells);
            evt.tstart(exclude) = [];
            evt.tend(exclude) = [];
            evt.nActiveCells(exclude) = [];

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
            
            % make "replay" position vector
%             spd_factor = 20;
%             path_len = 80; 
%             path_idx = pos_in.data > max(pos_in.data)-path_len;
%             pos_in2 = pos_in;
%             pos_in2.tvec = pos_in.tvec(path_idx)./spd_factor;
%             pos_in2.data = pos_in.data(path_idx);
            
%             cfg_model.nTrials = nSWR;
%             cfg_model.lockingParemeter = 0.7;
% %             model_run = pfmodel2(cfg_model,pfs,pos_in);
%             
% %             cfg_model.maxFiringRate = 5;
%             cfg_model.lockingParemeter = 0;
%             model_replay = pfmodel2(cfg_model,pfs,pos_in2);
% 
%             % do stuff to get "infield" data
% 
%             spikes = ConvertModel([],model_replay);
%             PlotSpikeRaster2([],spikes);
            
            % make place fields
            for i=1:length(field_order)
                pfs(i).xc = TC_temp(iT).field_loc(i);          
                pfs(i).maxFiringRate = data_in.IFFR_fieldmax(iT).(rn){1,iFD}(i);
            end
            
            curr_R = zeros(length(field_order),nSWR);
            spk_t = [];
            isi = [];
            
            for iTrial = nSWR:-1:1 %for each trial...
                pathlen = data_in.path_len(iT).(rn){1,iFD}(iTrial);
                if pathlen<1
                    continue
                end
                
                path_idx = pos_in.data > max(pos_in.data)-pathlen;
                
                tlen = data_in.SWR_size(iT).(rn){1,iFD}(iTrial);
                
                pos_in2 = pos_in;
                pos_in2.data = pos_in.data(path_idx);
                pos_in2.tvec = linspace(0,tlen,length(pos_in2.data));

                if isempty(pos_in2.data)
                    sprintf('hi')
                end
                
                cfg_model.nTrials = 1;
                cfg_model.lockingParemeter = 0;
                cfg_model.convertfile = 'spikes';
                curr_S = pfmodel2(cfg_model,pfs,pos_in2);

                for iC=1:length(curr_S.t) %for each cell...
                    curr_R(iC,iTrial) = length(curr_S.t{iC});
                    spk_t = vertcat(spk_t,curr_S.t{iC});
                    isi = vertcat(isi,diff(curr_S.t{iC}));

        %             % get "infield" data
        %             curr_tc = model.pf(iC).rates;
        %             [pks.loc,pks.full] = find_fields([],curr_tc);
        %             field_start = pks.full{:}(1);
        %             field_end = pks.full{:}(end);
        %             avg = mean(curr_tc(field_start:field_end));


                end %iterate cells
            end %iterate trials

            % get mean and variance across SWRs (for each unit)
            simdata.mean = vertcat(simdata.mean,mean(curr_R,2));
            simdata.var = vertcat(simdata.var,var(curr_R,0,2));

            % get length of each SWR event
            simdata.tlen = vertcat(simdata.tlen,max(spk_t)-min(spk_t));

            % get proportion of SWR each cell is active
            simdata.prop = vertcat(simdata.prop,sum(curr_R>0,2)/nSWR*100.0);

            % get mean spike count of place cells across SWR events
            curr_R(curr_R==0) = nan; %SWR events with no spikes are omitted from the average
            simdata.avgspk = vertcat(simdata.avgspk,nanmean(curr_R,2));

            % get ISIs
            simdata.ISI = vertcat(simdata.ISI,isi);
            
            clear curr_R curr_S
        end
    end
end

data_out = simdata;
profile viewer

end