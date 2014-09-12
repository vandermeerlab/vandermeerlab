%% load the data
clear all; pack
%cd('C:\Users\mvdm\Dropbox\teaching\CoSMo2014\R042-2013-08-18'); % isidro
%cd('D:\My_Documents\Dropbox\teaching\CoSMo2014\R042-2013-08-18'); % equinox
cd('D:\My_Documents\My Dropbox\teaching\CoSMo2014\R042-2013-08-18'); % athena

load(FindFile('*vt.mat'));
load(FindFile('*times.mat')); % this should go in metadata
load(FindFile('*CoordL.mat')); % this should go in metadata

cfg = [];
cfg.load_questionable_cells = 1;
S = LoadSpikes(cfg);

S = restrict(S,S.cfg.ExpKeys.TimeOnTrack,S.cfg.ExpKeys.TimeOffTrack);

%% set up spikes and position data for making tuning curves (encoding model)
ENC_S = restrict(S,run_start,run_end);
[ENC_S,keep] = removeEmptyCells(ENC_S);

S.t = S.t(keep); S.label = S.label(keep);

ENC_pos = restrict(pos,run_start,run_end);

%% obtain tuning curves; need to define bins first
SET_xmin = 80; SET_ymin = 0; % could go in metadata?
SET_xmax = 660; SET_ymax = 520;
SET_xBinSz = 10; SET_yBinSz = 10;

cfg = [];
cfg.binEdges{1} = SET_xmin:SET_xBinSz:SET_xmax;
cfg.binEdges{2} = SET_ymin:SET_yBinSz:SET_ymax;
cfg.smoothingKernel = gausskernel([4 4],2);

tc = TuningCurves(cfg,ENC_S,ENC_pos); 
%% plot example tc
iC = 7;
imagescnan(sq(tc.tc2D(iC,:,:))); axis off; colorbar;

%% plot all tc's
ppf = 25; % plots per figure
for iC = 1:length(ENC_S.t)
    nFigure = ceil(iC/ppf);
    figure(nFigure);
 
    subtightplot(5,5,iC-(nFigure-1)*ppf);
    h = imagescnan(sq(tc.tc2D(iC,:,:))); axis off;
    caxis([0 10]);
end

%% make Q-matrix
cfg = [];
cfg.dt = 0.25;
cfg.tvec_edges = run_start(1):cfg.dt:run_end(end);

Q = MakeQfromS(cfg,S);

% imagesc(Q.tvec,1:size(Q.data,1),Q.data) % plot
Q = restrict(Q,run_start,run_end); % TO FIX: RESTRICT KILLS USER FIELD

%% decode
cfg = [];
p = DecodeZ(cfg,Q,tc.tc);

%% plot
xBinned = interp1(ENC_pos.tvec,tc.pos_idx(:,1),Q.tvec);
yBinned = interp1(ENC_pos.tvec,tc.pos_idx(:,2),Q.tvec);

dec_err = nan(length(Q.tvec),1);

for iT = 1:length(Q.tvec)
    
    cla;
    
    toPlot = nan(tc.nBins{1},tc.nBins{2});
    toPlot(tc.good_idx) = p.data(iT,:);
    
    imagescnan(toPlot); axis xy; hold on; caxis([0 0.5]);
    axis off;
 
    hold on; plot(yBinned(iT),xBinned(iT),'ow','MarkerSize',15);
 
    % get x and y coordinates of bin with highest probability
    [~,idx] = max(toPlot(:));
    [x_map,y_map] = ind2sub(size(toPlot),idx);
    
    dec_err(iT) = sqrt((yBinned(iT)-y_map).^2+(xBinned(iT)-x_map).^2);
    
    plot(y_map,x_map,'g*','MarkerSize',5);
    
    h = title(sprintf('t %.2f, nCells %d, dist %.2f',Q.tvec(iT),Inf,dec_err(iT))); 

    drawnow; pause(0.1);
end


%% get distance (no plotting)
dec_err = nan(length(Q.tvec),1);

for iT = 1:length(Q.tvec);
    
    toPlot = nan(tc.nBins{1},tc.nBins{2});
    toPlot(tc.good_idx) = p.data(iT,:);
    
    [~,idx] = max(toPlot(:));
    [x_map,y_map] = ind2sub(size(toPlot),idx);
    
    dec_err(iT) = sqrt((yBinned(iT)-y_map).^2+(xBinned(iT)-x_map).^2);
    
end
plot(Q.tvec,dec_err,'.k');

%% plot better
% first, get trial id for each sample
trial_id = zeros(size(Q.tvec));
trial_idx = nearest_idx3(run_start,Q.tvec);
trial_id(trial_idx) = 1;
trial_id = cumsum(trial_id);

figure; set(gca,'FontSize',18);
boxplot(dec_err,trial_id);
xlabel('trial'); ylabel('decoding error (pixels)');

av_error = nanmean(dec_err);
title(sprintf('avg err %.2f',av_error));
