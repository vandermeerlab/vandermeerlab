%%
%SET_GitHub_root = 'C:\Users\mattm\Documents\GitHub'; % replace this with the location of your local GitHub root folder
SET_GitHub_root = 'C:\Users\mvdm\Documents\GitHub';
SET_data_fd = 'C:\data\R064\R064-2015-04-20'; % replace this with the location of your local data folder, download from https://drive.google.com/file/d/19lO0KQQjK9NQjr3xOz38bHwaG2yctr_q/view?usp=drive_link

restoredefaultpath;
addpath(genpath(cat(2,SET_GitHub_root,'\vandermeerlab\code-matlab\shared'))); % clone vandermeerlab repo at https://github.com/vandermeerlab/vandermeerlab

cfg_master = []; % holds config options for later in this script
%%
cd(SET_data_fd);
 
please = []; please.load_questionable_cells = 1;
S = LoadSpikes(please); % `please` variable overrides default LoadSpikes() options
 
pos = LoadPos([]); % empty input [] causes LoadPos() to use default options

LoadExpKeys; % annotation file containing some basic information about this data session
LoadMetadata; % loads experimenter-annotated file associated with each data session

%%
plot(getd(pos,'y'),getd(pos,'x'),'.','Color',[0.5 0.5 0.5],'MarkerSize',1); % note getd() gets data corresponding to a specific label (x and y here)
axis off; hold on;

iC = 7; % select cell 7 (out of 107 total)
spk_x = interp1(pos.tvec,getd(pos,'x'),S.t{iC},'linear');
spk_y = interp1(pos.tvec,getd(pos,'y'),S.t{iC},'linear');
 
h = plot(spk_y,spk_x,'.r');

%%
LoadMetadata; % loads experimenter-annotated file associated with each data session

% ENCoding variables: used for estimating tuning curves (encoding model)
ENC_S = restrict(S,metadata.taskvars.trial_iv); % trial_iv contains the start and end times of trials
ENC_pos = restrict(pos,metadata.taskvars.trial_iv);
 
% check for empties and remove
keep = ~cellfun(@isempty,ENC_S.t);
ENC_S = SelectTS([],ENC_S,keep);

% also set up DECoding variables for use later
DEC_S = SelectTS([],S,keep);

%%
figure;
cfg_tc = [];
cfg_tc.minOcc = 1;
cfg_tc.binEdges{1} = 50:10:610; cfg_tc.binEdges{2} = 0:10:480;
tc = TuningCurves(cfg_tc,ENC_S,ENC_pos);

pcolor(sq(tc.tc2D(iC,:,:))); shading flat; axis off; colorbar

%% try with smoothing
cfg_tc.smoothingKernel = gausskernel([4 4],2); % Gaussian kernel of 4x4 pixels, SD of 2 pixels (note this should sum to 1)
tc = TuningCurves(cfg_tc,ENC_S,ENC_pos);

pcolor(sq(tc.tc2D(iC,:,:))); shading flat; axis off; colorbar

%%
cfg_Q = [];
cfg_Q.binsize = 0.25;
cfg_Q.tvec_edges = metadata.taskvars.trial_iv.tstart(1):cfg_Q.binsize:metadata.taskvars.trial_iv.tend(end);
cfg_Q.tvec_centers = cfg_Q.tvec_edges(1:end-1)+cfg_Q.binsize/2;

Q = MakeQfromS(cfg_Q,DEC_S);
imagesc(Q.data);
%% decode
Q = restrict(Q,metadata.taskvars.trial_iv);

nBins = size(tc.tc, 2);
occUniform = repmat(1/nBins,[nBins 1]);

nActiveNeurons = sum(Q.data > 0);
 
len = length(Q.tvec);
p = nan(length(Q.tvec),nBins);
for iB = 1:nBins
    tempProd = nansum(log(repmat(tc.tc(:,iB),1,len).^Q.data));
    tempSum = exp(-cfg_Q.binsize*nansum(tc.tc(:,iB)',2));
    p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
end
 
p = p./repmat(sum(p,2),1,nBins); % renormalize to 1 total probability
p(nActiveNeurons < 1,:) = 0; % ignore bins with no activity

%%
xBinned = interp1(ENC_pos.tvec,tc.pos_idx(:,1),Q.tvec);
yBinned = interp1(ENC_pos.tvec,tc.pos_idx(:,2),Q.tvec);

%%
dec_err = nan(length(Q.tvec),1);
h = figure; set(h,'Position',[100 100 640 480]);

for iT = 1:length(Q.tvec)
    cla;
    toPlot = nan(tc.nBins{1},tc.nBins{2});
    toPlot(tc.good_idx) = p(iT,:);
 
    pcolor(toPlot); axis xy; hold on; caxis([0 0.5]);
    shading flat; axis off;
 
    hold on; plot(yBinned(iT),xBinned(iT),'ow','MarkerSize',15);
 
    % get x and y coordinates of MAP
    [~,idx] = max(toPlot(:));
    [x_map,y_map] = ind2sub(size(toPlot),idx);
 
    if nActiveNeurons(iT) > 0
        dec_err(iT) = sqrt((yBinned(iT)-y_map).^2+(xBinned(iT)-x_map).^2);
    end
 
    plot(y_map,x_map,'g*','MarkerSize',5);
 
    h = title(sprintf('t %.2f, nCells %d, dist %.2f',Q.tvec(iT),nActiveNeurons(iT),dec_err(iT))); 
    if nActiveNeurons(iT) == 0
        set(h,'Color',[1 0 0]);
    else
        set(h,'Color',[0 0 0]);
    end
    pause(0.1);
    % f(iT) = getframe(gcf); % store current frame to make video later
    drawnow;
end

%% optional- to make video
%fname = 'test.avi';
%vidObj = VideoWriter(fname,'MPEG-4'); vidObj.FrameRate = 10;
%open(vidObj);
%writeVideo(vidObj,f);
%close(vidObj);

%% error by space
cfg = [];
cfg.y_edges = cfg_tc.binEdges{2}; cfg.x_edges = cfg_tc.binEdges{1};
 
dec_err_tsd = tsd(Q.tvec,dec_err);
space_err = TSDbySpace(cfg,ENC_pos,dec_err_tsd);
 
figure;
pcolor(space_err); shading flat; axis off; colorbar; caxis([0 10]);

%% replay plotting
clear expCond;

expCond(1).label = 'left'; % this is a T-maze, we are interested in 'left' and 'right' trials
expCond(2).label = 'right'; % these are just labels we can make up here to keep track of which condition means what

expCond(1).t = metadata.taskvars.trial_iv_L; % previously stored trial start and end times for left trials
expCond(2).t = metadata.taskvars.trial_iv_R; 

expCond(1).coord = metadata.coord.coordL; % previously user input idealized linear track
expCond(2).coord = metadata.coord.coordR; % note, this is in units of "camera pixels", not cm

expCond(1).S = S;
expCond(2).S = S;

%% linearize paths (snap x,y position samples to nearest point on experimenter-drawn idealized track)
nCond = length(expCond);
for iCond = 1:nCond
   
    this_coord.coord = expCond(iCond).coord; this_coord.units = 'au'; this_coord.standardized = 0;
    expCond(iCond).linpos = LinearizePos([],pos,this_coord);
   
end

%% find intervals where rat is running
spd = getLinSpd([],pos); % get speed (in "camera pixels per second")

cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 10; 
run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above 10 pix/s


%% restrict (linearized) position data and spike data to desired intervals
for iCond = 1:nCond
   
    fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
end

%% get tuning curves
for iCond = 1:nCond
    
    cfg_tc = [];
    expCond(iCond).tc = TuningCurves(cfg_tc,expCond(iCond).S,expCond(iCond).linpos); 
    
end

%% detect place fields
for iCond = 1:nCond
    
    expCond(iCond).fields = DetectPlaceCells1D([],expCond(iCond).tc.tc);
    
end
    
%% load CSC with good SWRs
cfg = []; cfg.fc = ExpKeys.goodSWR(1);
lfp = LoadCSC(cfg); lfp = decimate_tsd([],lfp);

load(FindFile('*candidates*')); % SWR candidates previously saved, or can compute your own based on SWR_detection.m script

%% make rasterplot place cells for left and right separately, ordered by place field location
fh = figure('KeyPressFcn',@navigate);

for iCond = 1:nCond
    ax(iCond) = subplot(2,1,iCond);
    
    S_temp = SelectTS([],S,expCond(iCond).fields.template_idx); % template_idx contains ordered place cells
    
    cfg_mr = []; cfg_mr.openNewFig = 0; cfg_mr.lfp = lfp; cfg_mr.evt = evt; cfg_mr.lfpColor = 'k';
    MultiRaster(cfg_mr,S_temp); % see http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2016:week3short
    
    ylim([-25 40])
    xlim([6572 6574]); % zoom in on specific time interval of interest
end
linkaxes(ax,'x')

%% decode with moving window
addpath(genpath(cat(2,SET_GitHub_root,'\vandermeerlab\code-matlab\tasks\Alyssa_Tmaze')));

cfg = []; cfg.keepPosterior = 1; cfg.writeFiles = 0;
out = Generate_DecSeq(cfg);

%% visualize with decoding
ttl = {'left','right'}; 
xl = [6572 6574]; % R064-2015-04-20
for iCond = 1:nCond
    % raster
    ax(iCond) = subplot(2,2,iCond);
    out.expCond(iCond).fields = DetectPlaceCells1D([],out.expCond(iCond).tc.tc);
    S_temp = SelectTS([],out.expCond(iCond).decS,out.expCond(iCond).fields.template_idx);
    cfg_plot = [];
    %cfg_plot.lfp(1) = out.expCond(iCond).decode_map; 
    %cfg_plot.lfp(1) = lfp;
    %cfg_plot.evt = out.expCond(iCond).seq_iv;
    cfg_plot.openNewFig = 0;
    MultiRaster(cfg_plot,S_temp);
    %xlim(xl);
    set(gca,'YTick',[],'TickDir','out','XTick',[]);
    title(ttl{iCond}); xlabel(''); ylabel('');
    
    % decoded posterior
    ax(iCond+2) = subplot(2,2,iCond+2);
    imagesc(out.expCond(iCond).P2.tvec,1:111,out.expCond(iCond).P2.data);
    %set(gca,'YDir','normal','YTick',[],'XTick',xl(1):0.2:xl(2),'TickDir','out');
    set(gca,'YDir','normal','YTick',[],'TickDir','out');
    xlabel('time (s)');
    caxis([0 0.3]);
    xlim(xl);
end
linkaxes(ax,'x')