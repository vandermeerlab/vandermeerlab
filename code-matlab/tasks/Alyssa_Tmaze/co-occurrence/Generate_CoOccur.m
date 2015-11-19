function out = Generate_CoOccur(cfg_in)
% function out = Generate_CoOccur(cfg_in)
%
% run co-occurrence analysis from within session data folder
%
% assumes cfg_in.output_dir exists; image and data files are written into this
%
% CONFIGS:
%
% cfg_def = [];
% cfg_def.load_questionable_cells = 1;
% cfg_def.output_dir = 'files';
% cfg_def.output_file_prefix = []; % when writing files, prefix this
% cfg_def.whichEvents = 'all'; % 'prerecord', 'all', 'postrecord'
% cfg_def.dt = 0.1; % window size (s) for RUN
% cfg_def.win = 0.1; % window size (s) for SWR
% cfg_def.writeFiles = 1; % write image (if selected) and data files
% cfg_def.plotOutput = 0; % generate images
% cfg_def.nShuffles = 250; % number of shuffles for co-occurrence z-scoring
% cfg_def.removeTheta = 1; % remove events with high theta
% cfg_def.removeRunning = 1; % remove events with high running speed
% cfg_def.doMovingWindow = 0; % perform moving window analysis
%
% MvdM 2015

cfg_def = [];
cfg_def.load_questionable_cells = 1;
cfg_def.output_dir = 'files';
cfg_def.output_file_prefix = []; % when writing files, prefix this
cfg_def.whichEvents = 'all'; % 'prerecord', 'all', 'postrecord', 'taskrest'
cfg_def.dt = 0.1; % window size (s) for RUN
cfg_def.win = 0.1; % window size (s) for SWR
cfg_def.writeFiles = 1;
cfg_def.plotOutput = 0;
cfg_def.nShuffles = 250;
cfg_def.removeTheta = 1;
cfg_def.removeRunning = 1;
cfg_def.doMovingWindow = 0;
cfg_def.outputPConly = 0; % send me the place cells and tuning curves only

cfg = ProcessConfig2(cfg_def,cfg_in);

LoadMetadata
LoadExpKeys

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

% set up filenames and paths
this_fd = pwd;
output_fd = cat(2,pwd,'\',cfg.output_dir);
base_fn = cat(2,cfg.output_file_prefix,S_orig.cfg.SessionID);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% linearization block %%
%%%%%%%%%%%%%%%%%%%%%%%%%
pos = LoadPos([]);

% get similar trial numbers
if strcmp(S_orig.cfg.SessionID(1:4),'R042')
    metadata = TrimTrialTimes([],metadata); % R042 only!!
end
%[L_trl,R_trl] = GetMatchedTrials([],metadata,ExpKeys);
[L_trl,R_trl] = GetMatchedTrials_old([],metadata);

% restrict spiketrains to trial times only; note that some spiketrains may
% be empty after this step
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

%ENC is for "encoding" ie the tuning curves, that describe the
%relationship between our observable variable (position) and firing rate.
%(the answer to a question like, "How is position (en)coded by the
%neurons whose activity we are recording?")
ENC_data(1).trial_type = 'left';
ENC_data(1).Coord = CoordLrs;
ENC_data(1).pos = posL_binned;
ENC_data(1).S = S_left; % all spikes from all units during left trials

ENC_data(2).trial_type = 'right';
ENC_data(2).Coord = CoordRrs;
ENC_data(2).pos = posR_binned;
ENC_data(2).S = S_right; % all spikes from all units during right trials

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

% write fig
if cfg.writeFiles && cfg.plotOutput, cd(output_fd);
    out_fn = cat(2,base_fn,'-CoordCheck.png');
    print(gcf,'-r300','-dpng',out_fn);
    close(gcf);
    cd(this_fd);
end

clear CoordL CoordLrs CoordR CoordRrs S_left S_right posL posL_binned posR posR_binned nBins trial_ID t

%% apply speed filter to encoding data
spd = getLinSpd([],pos);
cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 10; run_iv = TSDtoIV(cfg_spd,spd);

% now the spiketrains must contain only spikes that were emitted when the rat
% was traveling faster than 10 pixels/s:
ENC_data(1).S = restrict(ENC_data(1).S,run_iv);
ENC_data(2).S = restrict(ENC_data(2).S,run_iv);

% keep this for later display
S_run = restrict(S_orig,metadata.taskvars.trial_iv);

%% Make tuning curves for left or right trajectory
cfg_tc = [];
cfg_tc.binSize = 1; % now in units of pos...

TC_left = MakeTC(cfg_tc,ENC_data(1).S,ENC_data(1).pos);
TC_right = MakeTC(cfg_tc,ENC_data(2).S,ENC_data(2).pos);

cfg_tc = []; cfg_tc.order = 2; cfg_tc.binsize = binSize; cfg_tc.cp = cpL;% fields (1 is peaks)

if cfg.plotOutput, TCPlot(cfg_tc,TC_left); end

% write fig -- could be helper function
if cfg.writeFiles && cfg.plotOutput, cd(output_fd);
    out_fn = cat(2,base_fn,'-TC_left.png');
    print(gcf,'-r300','-dpng',out_fn);
    close(gcf);
    cd(this_fd);
end

cfg_tc.cp = cpR;
if cfg.plotOutput, TCPlot(cfg_tc,TC_right); end
% write fig -- could be helper function
if cfg.writeFiles && cfg.plotOutput, cd(output_fd);
    out_fn = cat(2,base_fn,'-TC_right.png');
    print(gcf,'-r300','-dpng',out_fn);
    close(gcf);
    cd(this_fd);
end

% NOTE this takes a lot of Java heap space..
if cfg.plotOutput, fh = AnotherTCplotTwin([],S_run,TC_left,TC_right,pos); else fh = []; end
% write figs
if cfg.writeFiles && cfg.plotOutput, cd(output_fd);
    for iFH = 1:length(fh)
        out_fn = cat(2,base_fn,'-fields_',num2str(iFH),'.png');
        set(fh(iFH), 'InvertHardCopy', 'off');
        print(fh(iFH),'-r300','-dpng',out_fn);
        close(fh(iFH));
    end
    cd(this_fd);
end

%% order data -- contents of S are now fields beyond the choice point
FieldOrder_left = TC_left.field_template_idx(TC_left.field_loc > cpL); % this keeps all cells with field after the choice point
FieldOrder_right = TC_right.field_template_idx(TC_right.field_loc > cpR); % ^^^

% remove cells with fields in both
[~,ia,ib] = intersect(FieldOrder_left,FieldOrder_right);
FieldOrder_left(ia) = []; FieldOrder_right(ib) = [];

% make a new S, keeping left-only cells
S_pc_left = OrderSelectS([],S_orig,FieldOrder_left);
%S_pc_left.t = S_pc_left.t(FieldOrder_left); S_pc_left.label = S_pc_left.label(FieldOrder_left); S_pc_left.usr.data = S_pc_left.usr.data(FieldOrder_left);

% make a new S, keeping right-only cells
S_pc_right = OrderSelectS([],S_orig,FieldOrder_right);
%S_pc_right.t = S_pc_right.t(FieldOrder_right); S_pc_right.label = S_pc_right.label(FieldOrder_right); S_pc_right.usr.data = S_pc_right.usr.data(FieldOrder_right);

if cfg.outputPConly
    % output data pertaining to the place cells we keep for analysis
    
    keepLeft = TC_left.field_template_idx(TC_left.field_loc > cpL);
    keepRight = TC_right.field_template_idx(TC_right.field_loc > cpR);
    
    %spiketrains and corresponding .t file name
    
    armPC.L.S = OrderSelectS([],S_orig,TC_left.template_idx);
    %tuning curves
    armPC.L.tc = TC_left; %%%% was here
    armPC.L.unique_idx = keepLeft;
    
    armPC.R.S = OrderSelectS([],S_orig,TC_right.template_idx);
    
    armPC.R.tc = TC_right;
    armPC.R.unique_idx = keepRight;
    cd(output_fd)
    save([base_fn,'-armPC'],'armPC')
    cd(this_fd)
    %return
end
if ~cfg.outputPConly
    %% make Q-matrix for in-field activity
    S_lr = []; % make a new S containing left-only and right-only cells
    S_lr.t = cat(2,S_pc_left.t,S_pc_right.t); nLeft = length(S_pc_left.t); nRight = length(S_pc_right.t);
    S_lr = concatenateTS(S_pc_left,S_pc_right);
    
    %A Q-matrix is organized such that each row corresponds to a cell, and each
    %column groups the spikes into time bins. Thus, each column contains information
    %about which cells were active together in a given time bin.
    cfg_Q = [];
    cfg_Q.dt = cfg.dt;
    cfg_Q.tvec_edges = ExpKeys.prerecord(1):cfg_Q.dt:ExpKeys.postrecord(end);
    
    Q_lr = MakeQfromS(cfg_Q,S_lr);
    
    % restrict to trials/running
    both_trl = iv;
    both_trl.tstart = cat(1,L_trl.tstart,R_trl.tstart);
    both_trl.tend = cat(1,L_trl.tend,R_trl.tend);
    Q_lrR = restrict(Q_lr,both_trl);
    
    spd = getLinSpd([],pos);
    cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 20; run_iv = TSDtoIV(cfg_spd,spd);
    Q_lrR = restrict(Q_lrR,run_iv);
    
    % restrict to maze pieces beyond choice point...
    cfg_l = []; cfg_l.method = 'raw'; cfg_l.threshold = cpL; cfg_l.target = 'z'; run_left = TSDtoIV(cfg_l,ENC_data(1).pos);
    cfg_r = []; cfg_r.method = 'raw'; cfg_r.threshold = cpR; cfg_r.target = 'z'; run_right = TSDtoIV(cfg_r,ENC_data(2).pos);
    run_iv = UnionIV([],run_left,run_right);
    Q_lrR = restrict(Q_lrR,run_iv);
    
    %% z-scored co-occurrence
    cx = [-20 20];
    
    cfg_cc = []; cfg_cc.nShuffle = cfg.nShuffles; cfg_cc.num_tt = S_lr.usr.tt_num;
    [c_crossed,mask_both] = CoOccurQ(cfg_cc,Q_lrR);
    
    if cfg.plotOutput
        figure;
        imagesc(c_crossed.p4); colorbar;
        hold on; plot(xlim,[nLeft+0.5 nLeft+0.5],'k--'); plot([nLeft+0.5 nLeft+0.5],ylim,'k--');
        tl = [FieldOrder_left FieldOrder_right];
        set(gca,'XTick',1:length(tl),'XTickLabel',tl,'YTick',1:length(tl),'YTickLabel',tl);
        title(sprintf('RUN - coOcc z (dt = %.3f)',cfg.dt)); caxis(cx);
    end
    % write fig -- could be helper function
    if cfg.writeFiles && cfg.plotOutput, cd(output_fd);
        out_fn = cat(2,base_fn,'-RUN_coOccurZ.png');
        print(gcf,'-r300','-dpng',out_fn);
        close(gcf);
        cd(this_fd)
    end
    
    % some stats
    all_left = c_crossed.p4(1:nLeft,1:nLeft);
    all_right = c_crossed.p4(nLeft+1:end,nLeft+1:end);
    all_crossed = c_crossed.p4(1:nLeft,nLeft+1:end);
    
    fprintf('\nRUN left mean z %.2f +/- %.2f, median %.2f\n',nanmean(all_left(:)),nanstd(all_left(:)),nanmedian(all_left(:)));
    fprintf('RUN right mean z %.2f +/- %.2f, median %.2f\n',nanmean(all_right(:)),nanstd(all_right(:)),nanmedian(all_right(:)));
    fprintf('RUN crossed mean z %.2f +/- %.2f, median %.2f\n',nanmean(all_crossed(:)),nanstd(all_crossed(:)),nanmedian(all_crossed(:)));
    
    % keep
    out.run.full_cc = c_crossed;
    out.run.left = all_left;
    out.run.right = all_right;
    out.run.crossed = all_crossed;
    out.nLeft = nLeft; out.nRight = nRight;
    out.TCleft = TC_left; out.TCright = TC_right;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% get SWR candidates %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %load(FindFile('*candidates.mat'));
    LoadCandidates
    if isfield(evt,'data'), evt = rmfield(evt,'data'); end
    
    fprintf('nEvents loaded: %d\n',length(evt.tstart));
    
    switch cfg.whichEvents
        case 'all'
            % do nothing
        case 'prerecord'
            evt = restrict(evt,ExpKeys.prerecord(1),ExpKeys.prerecord(2));
        case 'postrecord'
            evt = restrict(evt,ExpKeys.postrecord(1),ExpKeys.postrecord(2));
        case 'taskrest'
            evt = restrict(evt,metadata.taskvars.trial_iv.tstart,metadata.taskvars.trial_iv.tend);
        case 'task'
            evt = restrict(evt,ExpKeys.task(1),ExpKeys.task(2));
        otherwise
            error('Unrecognized cfg.whichEvents')
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% co-coccurrence per SWR %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % make a Q-matrix
    %A Q-matrix is organized such that each row corresponds to a cell, and each
    %column groups the spikes into time bins. Thus, each column contains information
    %about which cells were active together in a given time bin.
    % the time bins here are +/- 50 ms of the center time of the candidate
    % events; thus, each time bin is 100 ms wide
    
    cfg_Q = [];
    
    nEvt = length(evt.tend);
    Q_leftSWR = []; Q_rightSWR = [];
    
    if ~isempty(S_pc_left.t)
        for iEvt = nEvt:-1:1
            
            %cfg_Q.tvec_edges = [evt.tstart(iEvt) evt.tend(iEvt)];
            center_t = mean([evt.tstart(iEvt) evt.tend(iEvt)]);
            cfg_Q.tvec_edges = [center_t-cfg.win/2 center_t+cfg.win/2];
            cfg_Q.dt = diff(cfg_Q.tvec_edges);
            
            qtemp = MakeQfromS(cfg_Q,S_pc_left);
            Q_leftSWR(:,iEvt) = qtemp.data;
            
        end
    end
    
    if ~isempty(S_pc_right.t)
        for iEvt = nEvt:-1:1
            
            %cfg_Q.tvec_edges = [evt.tstart(iEvt) evt.tend(iEvt)];
            center_t = mean([evt.tstart(iEvt) evt.tend(iEvt)]);
            cfg_Q.tvec_edges = [center_t-cfg.win/2 center_t+cfg.win/2];
            cfg_Q.dt = diff(cfg_Q.tvec_edges);
            
            qtemp = MakeQfromS(cfg_Q,S_pc_right);
            Q_rightSWR(:,iEvt) = qtemp.data;
            assignin('base','Q',Q_rightSWR)
            
        end
    end
    
    Q_leftSWRq.data = Q_leftSWR;
    Q_rightSWRq.data = Q_rightSWR;
    
    %% z-scored
    cfg_cc = []; cfg_cc.num_tt = cat(2,S_pc_left.usr.tt_num,S_pc_right.usr.tt_num);
    cfg_cc.nShuffle = cfg.nShuffles;
    
    Q_leftright_SWRq = [];
    Q_leftright_SWRq.data = cat(1,Q_leftSWRq.data,Q_rightSWRq.data);
    
    % "crossed" is the term for pairs with one cell on the left
    %arm and another on the right arm -- so those we would expect to be less
    %co-active (at least during behavior) than pairs of cells on the same arm.
    c_crossed = CoOccurQ(cfg_cc,Q_leftright_SWRq);
    
    if cfg.plotOutput
        figure;
        cx = [-5 5];
        imagesc(c_crossed.p4); colorbar;
        hold on; plot(xlim,[nLeft+0.5 nLeft+0.5],'k--'); plot([nLeft+0.5 nLeft+0.5],ylim,'k--');
        tl = [FieldOrder_left FieldOrder_right];
        set(gca,'XTick',1:length(tl),'XTickLabel',tl,'YTick',1:length(tl),'YTickLabel',tl);
        title(sprintf('SWR_z - coOcc z (dt = %.3f)',cfg.dt)); caxis(cx);
    end
    % write fig -- could be helper function
    if cfg.writeFiles && cfg.plotOutput, cd(output_fd);
        out_fn = cat(2,base_fn,'-coOccur_SWR_z.png');
        print(gcf,'-r300','-dpng',out_fn);
        close(gcf);
        cd(this_fd);
    end
    
    % some stats
    all_left = c_crossed.p4(1:nLeft,1:nLeft); % nLeft is number of left cells. this is keeping stat data made from the upper left quadrant of the Q-matrix, which is the left-only cells
    %since left and then right data are concatenated, the first nLeft entries are for left cells, and the rest are for right cells
    all_right = c_crossed.p4(nLeft+1:end,nLeft+1:end); % this is keeping the stat data made from lower right quadrant of the Q-matrix, which is the right-only cells
    all_crossed = c_crossed.p4(1:nLeft,nLeft+1:end); % this is keeping stat data from the ?upper right?, which is for
    % left-right pairs. the remaining quadrant should be indentical to this one because it's a symmetrical matrix?
    
    fprintf('\nSWRz left mean z %.2f +/- %.2f, median %.2f\n',nanmean(all_left(:)),nanstd(all_left(:)),nanmedian(all_left(:)));
    fprintf('SWRz right mean z %.2f +/- %.2f, median %.2f\n',nanmean(all_right(:)),nanstd(all_right(:)),nanmedian(all_right(:)));
    fprintf('SWRz crossed mean z %.2f +/- %.2f, median %.2f\n',nanmean(all_crossed(:)),nanstd(all_crossed(:)),nanmedian(all_crossed(:)));
    
    % keep
    out.SWRz.full_cc = c_crossed;
    out.SWRz.left = all_left;
    out.SWRz.right = all_right;
    out.SWRz.crossed = all_crossed;
    
    %% moving window version
    if cfg.doMovingWindow
        winsize = 100;
        dw = 25;
        win_start = 1:dw:nEvt-winsize;
        win_end = win_start+winsize-1;
        win_centers = win_start+winsize/2;
        
        %clear l_mean r_mean c_mean l_sd r_sd c_sd;
        for iW = length(win_start):-1:1
            
            fprintf('Moving window %d/%d...\n',iW,length(win_start));
            
            this_Q = Q_leftright_SWRq;
            this_Q.data = this_Q.data(:,win_start(iW):win_end(iW));
            
            c_crossed = CoOccurQ(cfg_cc,this_Q);
            
            all_left = c_crossed.p4(1:nLeft,1:nLeft);
            all_right = c_crossed.p4(nLeft+1:end,nLeft+1:end);
            all_crossed = c_crossed.p4(1:nLeft,nLeft+1:end);
            
            out.l_mean(iW) = nanmean(all_left(:)); out.l_sd(iW) = nanstd(all_left(:));
            out.r_mean(iW) = nanmean(all_right(:)); out.r_sd(iW) = nanstd(all_right(:));
            out.c_mean(iW) = nanmean(all_crossed(:)); out.c_sd(iW) = nanstd(all_crossed(:));
            
            out.tvec_centerEvent = evt.tstart(win_centers);
            out.tvec_firstEvent = evt.tstart(win_start);
            out.tvec_lastEvent = evt.tstart(win_end);
            
        end
    end
    
    out.cfg = cfg;
    
    %% write data -- could be helper function
    if cfg.writeFiles
        cd(output_fd);
        out_fn = cat(2,base_fn,'-coOccur_data.mat');
        save(out_fn,'out');
        cd(this_fd)
    end
end
cd(this_fd)
