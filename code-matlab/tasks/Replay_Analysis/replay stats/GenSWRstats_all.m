function [data_out] = GenSWRstats_all(cfg_in)

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

iFD = 1; iR=1;
nCells_out = [];
nSWR_out = [];

for iR = 1:length(cfg.rats)
    rn = cfg.rats{iR};
    R(2).(rn) = [];
    SWR_size(2).(rn) = [];

    IFFR(2).(rn) = [];
    IFFR_field(2).(rn) = [];
    IFFR_fieldmax(2).(rn) = [];
    
    ISI_SWR(2).(rn) = [];
    ISI_IF(2).(rn) = [];
    
    TC.(rn) = [];
    ENC.(rn) = [];
    
    path_length(2).(rn) = [];
    
    cfg_temp.rats = cfg.rats(iR);
    fd = sort(getTmazeDataPath(cfg_temp)); % get all session directories  
    
    %% run for each session
    for iFD = 1:length(fd)
        dirname = fd{iFD};
        cd(dirname);

        LoadExpKeys
        LoadMetadata % for freqs

        %load position
        pos = LoadPos([]);

        %load good SWR csc
        cfg_temp.fc = ExpKeys.goodSWR(1);
        csc = LoadCSC(cfg_temp); 

        %load candidate events
        load(FindFile('*candidates.mat'));
        if strcmp(cfg.restrict,'pre');
            evt = restrict2(evt,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
        elseif strcmp(cfg.restrict,'post');
            evt = restrict2(evt,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
        elseif strcmp(cfg.restrict,'task');
            evt = restrict2(evt,ExpKeys.task(1),ExpKeys.task(2));
        end
        
        % Load spikes
        cfg_S = [];
        cfg_S.useClustersFile = 0;
        cfg_S.load_questionable_cells = 0;
        cfg_S.min_cluster_quality = 3;
        cfg_S.getTTnumbers = 1;
        S_orig = LoadSpikes(cfg_S); clear cfg_S
        
        %% Get place cells
        % get similar trial numbers
        if strcmp(S_orig.cfg.SessionID(1:4),'R042')
            metadata = TrimTrialTimes([],metadata); % R042 only!!
        end
        [L_trl,R_trl] = GetMatchedTrials([],metadata,ExpKeys);

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
            plot(ENC_data(1).pos.tvec,getd(ENC_data(1).pos,'z'))
            title('Left-linearized pos');
            subplot(224);
            plot(ENC_data(2).pos.tvec,getd(ENC_data(2).pos,'z'))
            title('Right-linearized pos');    
            pause(); close;
        end
        clear CoordL CoordLrs CoordR CoordRrs S_left S_right posL posL_binned posR posR_binned nBins trial_ID t

        % apply speed filter to encoding data
        spd = getLinSpd([],pos);
        cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 10; run_iv = TSDtoIV(cfg_spd,spd);

        ENC_data(1).S = restrict2(ENC_data(1).S,run_iv);
        ENC_data(2).S = restrict2(ENC_data(2).S,run_iv);
        ENC.(rn){iFD} = ENC_data;

        %% Make tuning curves for left or right trajectories
        clear TC_temp %in case you are re-running
        cfg_tc = [];
        cfg_tc.binSize = 1;
        TC_temp(1) = MakeTC(cfg_tc,ENC_data(1).S,ENC_data(1).pos);
        TC_temp(2) = MakeTC(cfg_tc,ENC_data(2).S,ENC_data(2).pos);
        TC.(rn){iFD} = TC_temp;
        
        cfg_tc = []; 
        cfg_tc.order = 2; 
        cfg_tc.binsize = binSize; 

        cfg_tc.cp = cpL;% fields (1 is peaks)
        if cfg.plotOutput, TCPlot(cfg_tc,TC_temp(1)); end

        cfg_tc.cp = cpR;
        if cfg.plotOutput, TCPlot(cfg_tc,TC_temp(2)); pause(); close all; end


        %% Generate stats for each arm
        run_trl = [L_trl R_trl];
        for iT = 1:2
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
                    
            % Generate Q-matrix
            cfg_temp = [];
            cfg_temp.dt = 0.005; %0.005 for NAU, 0.0005 for csc sample rate
            Q = MakeQfromS(cfg_temp,S_pc(iT));
            
            if cfg.NAU
                % limit candidates by NAU of place cells
                minCells = 4;
                activeCellsIV = AddNActiveCellsIV(cfg_temp,iv_in,S_pc(iT));
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


            % Generate R-matrix (nCell x nSWR matrix containing spike count in each SWR event interval)
            R_temp = zeros(nCells,nSWR);
            SWR_size_temp = zeros(nSWR,1);
            ISI_SWR_temp = [];
            path_length_temp = [];
            
            % for each SWR event
            for iIV = 1:nSWR
                % get spike counts for each cell for this event
                Qr = restrict2(Q,evt.tstart(iIV),evt.tend(iIV));
                R_temp(:,iIV) = sum(Qr.data,2);
     
                % for each place cell get ISI distribution
                S_SWR = restrict2(S_pc(iT),evt.tstart(iIV),evt.tend(iIV)); %restrict2
                for iC=1:nCells %for each place cell...
                    ISI_SWR_temp = vertcat(ISI_SWR_temp,diff(S_SWR.t{iC}));
                end
                                
                % for each SWR event, get how long the event is
                SWR_size_temp(iIV) = evt.tend(iIV)-evt.tstart(iIV);
                
                % get virtual path length of each event (% max dist between first and last peak)
                path_inds = find(R_temp(:,iIV));
                if isempty(path_inds);
                    path_length_temp = vertcat(path_length_temp,nan); 
                    continue
                end
%                 assert(~isempty(path_inds),'SWR event has no spikes')
                peak_loc_temp = [];
                for iP = 1:length(path_inds) % for each cell find peak
                    peak_loc_temp = horzcat(peak_loc_temp,TC_temp(iT).peak_loc{path_inds(iP)});
                end
                path_start = min(peak_loc_temp); 
                path_end = max(peak_loc_temp);
                path_length_temp = vertcat(path_length_temp,(path_end-path_start)); 

            end
            
            SWR_size(iT).(rn){iFD} = SWR_size_temp;
            R(iT).(rn){iFD} = R_temp;
            ISI_SWR(iT).(rn){iFD} = ISI_SWR_temp;
            path_length(iT).(rn){iFD} = path_length_temp;

            %% Get in-field firing rates only
            field_loc = TC_temp(iT).field_loc; 
            pos_mat = getd(ENC_data(iT).pos,'z');
            ISI_IF_temp = [];
            IFFR_field_temp = [];
            IFFR_fieldmax_temp = [];
            IFFR_temp = [];
            
            for iC=1:nCells; %for each place cell...
                curr_tc = TC_temp(iT).tc(:,field_order(iC));
                [pks.loc, pks.full] = find_fields(cfg_tc,curr_tc); %get field values
                
                PCFR_TEMP = []; %mean firing rate for place cell across peaks
                for curr_peak=1:length(pks.loc) %for each peak...
                    % get place field boundaries
                    field_start = pks.full{curr_peak}(1);
                    field_end = pks.full{curr_peak}(end);
                    
                    % get infield times
                    cfg_infield = [];
                    cfg_infield.method = 'raw'; cfg_infield.threshold = [field_start-1 field_end+1]; 
                    cfg_infield.dcn = 'range'; cfg_infield.target = 'z';
                    infield_iv = TSDtoIV(cfg_infield,ENC_data(iT).pos);
                    
                    % get infield ISIs
                    for iP=1:length(infield_iv.tstart) %for each pass of the field...
                        S_pass = restrict2(S_pc(iT),infield_iv.tstart(iP),infield_iv.tend(iP)); %restrict2
                        ISI_IF_temp = vertcat(ISI_IF_temp,diff(S_pass.t{iC}));
                    end
                    
                    % get mean and max firing rate during each field separately
                    IFFR_curr_field = mean(curr_tc(field_start:field_end));
                    IFFR_max_field = max(curr_tc(field_start:field_end));
                    IFFR_field_temp = vertcat(IFFR_field_temp,IFFR_curr_field);
                    IFFR_fieldmax_temp = vertcat(IFFR_fieldmax_temp,IFFR_max_field);

                    PCFR_TEMP = horzcat(PCFR_TEMP,IFFR_curr_field);
                    
                end %iterate peaks
                
                IFFR_temp = vertcat(IFFR_temp,mean(PCFR_TEMP));
            end %iterate place cells
            
            ISI_IF(iT).(rn){iFD} = ISI_IF_temp;
            IFFR(iT).(rn){iFD} = IFFR_temp;
            IFFR_field(iT).(rn){iFD} = IFFR_field_temp;
            IFFR_fieldmax(iT).(rn){iFD} = IFFR_fieldmax_temp;
            
            nCells_out = vertcat(nCells_out,nCells);
            nSWR_out = vertcat(nSWR_out,nSWR);
        end %iterate left vs right
        
    end %iterate sessions
    
end %iterate rats

data_out.TC = TC;
data_out.ENC_data = ENC;
data_out.R = R;
data_out.SWR_size = SWR_size;
data_out.IFFR = IFFR;
data_out.IFFR_field = IFFR_field;
data_out.IFFR_fieldmax = IFFR_fieldmax;
data_out.ISI_SWR = ISI_SWR;
data_out.ISI_IF = ISI_IF;
data_out.path_len = path_length;
data_out.nCells = nCells_out;
data_out.nSWR = nSWR_out;

cd(originalFolder)

end