%%
cd('D:\data\R066\R066-2014-11-27');
nTargets = 2;

%% load raw data
cfg = [];
cfg.targets = 1:nTargets;
cfg.mode = 'raw';

[ x_avg,y_avg,t_avg ] = xy_targets_chunks(cfg);
t_avg = t_avg * 10^-6; % convert to s

if ~exist('temp_pos.mat','file'); % speed up processing later
    save('temp_pos.mat','x_avg','y_avg','t_avg');
end

%% plot each target individually (raw)
figure(1)

subplot(221);
cols = 'rgbcm';
for iNT = nTargets:-1:1
    plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end
title('x');

subplot(222);
cols = 'rgbcm';
for iNT = nTargets:-1:1
    plot(t_avg,y_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end
title('y');

subplot(223);
hist(x_avg,1:720);
title('x');

subplot(224);
hist(y_avg,1:640);
title('y');

%% find intervals with LEDs on

cfg = []; 
cfg.eventList = {'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0001).', ...
    'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0002).', ...
    'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0000).'};

evt = LoadEvents(cfg);

led0 = iv(evt.t{1},evt.t{3}(nearest_idx3(evt.t{1},evt.t{3},1)));
led1 = iv(evt.t{2},evt.t{3}(nearest_idx3(evt.t{2},evt.t{3},1)));

%% check for nonnegative intervals
figure;
subplot(121)
hist(led0.tend-led0.tstart,100)
subplot(122)
hist(led1.tend-led1.tstart,100)

%% expand IVs -- this is apparently needed to catch stray position samples which fall outside the LED event window
cfg = []; cfg.d = [-1/60 1/60];
led0 = ResizeIV(cfg,led0);
led1 = ResizeIV(cfg,led1);

%% restrict (slow, need to optimize TSD_getidx)
xtsd = tsd(t_avg,x_avg);
led0_idx = TSD_getidx(xtsd,led0);
led1_idx = TSD_getidx(xtsd,led1);

%% plot each target individually (raw) plus points to be removed (x)
figure;
cols = 'rgbcm';
for iNT = nTargets:-1:1
    plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end
plot(t_avg(led0_idx),x_avg(led0_idx,1),'o','Color',cols(1));
plot(t_avg(led1_idx),x_avg(led1_idx,1),'o','Color',cols(1));

%% plot each target individually (raw) plus points to be removed (y)
figure;
cols = 'rgbcm';
for iNT = nTargets:-1:1
    plot(t_avg,y_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end
plot(t_avg(led0_idx),y_avg(led0_idx,1),'o','Color',cols(1));
plot(t_avg(led1_idx),y_avg(led1_idx,1),'o','Color',cols(1));

%% define LED centroids (take from plots)
% looks like: x0 = 521 +/- 5, x1 = 623 +/- 5; y0 = 473 +/- 5, y1 = 57 +/- 5
led0_xy_centroid = [521 473];
led1_xy_centroid = [623 57];

%% get distance of each sample to LEDs
DistFun = @(xd,yd,ctr) sqrt((xd-ctr(1)).^2 + (yd-ctr(2)).^2); % euclidian distance is only accurate if pixels are square (which they are not)

for iT = nTargets:-1:1
    led0_dist(:,iT) = DistFun(x_avg(:,iT),y_avg(:,iT),led0_xy_centroid);
    led1_dist(:,iT) = DistFun(x_avg(:,iT),y_avg(:,iT),led1_xy_centroid);
end

%% set position sample to nan if LED is on AND position is within radius
thr = 10;

% binarize LED on/off vectors for easier comparison
led0_bin = zeros(size(led0_dist(:,1))); led0_bin(led0_idx) = 1;
led1_bin = zeros(size(led1_dist(:,1))); led1_bin(led1_idx) = 1;

for iT = nTargets:-1:1
    
    nan_idx(:,iT) = (led0_bin & led0_dist(:,iT) < thr) | (led1_bin & led1_dist(:,iT) < thr);
    
    %x_avg(nan_idx(:,iT),iT) = NaN; y_avg(nan_idx(:,iT),iT) = NaN; 

end

%% plot result for targets 1-2
figure;
cols = 'rgbcm';
for iNT = 2:-1:1
    plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end
plot(t_avg(nan_idx(:,1)),x_avg(nan_idx(:,1),1),'o','Color',cols(1));
plot(t_avg(nan_idx(:,2)),x_avg(nan_idx(:,2),2),'o','Color',cols(2));

%% remove
for iNT = nTargets:-1:1
    
    x_avg(nan_idx(:,iNT),iNT) = NaN; y_avg(nan_idx(:,iNT),iNT) = NaN; 

end

%% plot result for targets 1-2
figure;
cols = 'rgbcm';
for iNT = 2:-1:1
    plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end

%% remove samples which jump + are in LED locations

jump_thr = 100;
all_removed = 0; iI = 1;

while ~all_removed
    
    fprintf('Removing jumps, iteration %d...\n',iI);
    
    clear nan_idx dist_diff;
    for iNT = 2:-1:1
        % first get distances between successive samples (diffs)
        dist_diff(:,iNT) = cat(1,0,sqrt(diffskipnan(x_avg(:,iNT)).^2 + diffskipnan(y_avg(:,iNT)).^2));
        
        nan_idx(:,iNT) = (dist_diff(:,iNT) > jump_thr) & ((led0_dist(:,iNT) < thr) | (led1_dist(:,iNT) < thr));
        
        x_avg(nan_idx(:,iNT),iNT) = NaN; y_avg(nan_idx(:,iNT),iNT) = NaN;
        
    end
    
    if all(nan_idx(:) == 0)
        all_removed = 1;
    end
    iI = iI + 1;
end

%% plot final result for targets 1-2
figure;
cols = 'rgbcm';
ax(1) = subplot(211);
for iNT = 2:-1:1
    plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end
plot(t_avg(nan_idx(:,1)),x_avg(nan_idx(:,1),1),'o','Color',cols(1));
plot(t_avg(nan_idx(:,2)),x_avg(nan_idx(:,2),2),'o','Color',cols(2));

ax(2) = subplot(212);
for iNT = 2:-1:1
    plot(t_avg,y_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end
plot(t_avg(nan_idx(:,1)),y_avg(nan_idx(:,1),1),'o','Color',cols(1));
plot(t_avg(nan_idx(:,2)),y_avg(nan_idx(:,2),2),'o','Color',cols(2));

linkaxes(ax,'x');

%% ad hoc method for reasonable tracking:
% (1) use red target if available
% (2) if not available, use green target if available
% (3) use a median filter to get rid of sporadic jumpy points (tried various
% Kalmans, but can't seem to get the right movement/measurement models for
% good performance)
% (4) purposefully did not try to interpolate missing data (can be done
% separately with e.g. interp1 spline, or Kalman

%% (1)
% take red target as starting point
x = x_avg(:,1); y = y_avg(:,1);

nan_idx = isnan(x);
x(nan_idx) = x_avg(nan_idx,2);
y(nan_idx) = y_avg(nan_idx,2);

nan_idx = isnan(x);

% filter
x = medfilt1m(x,5);
y = medfilt1m(y,5);

% put NaNs back
x(nan_idx) = nan; y(nan_idx) = nan;

axes(ax(1));
plot(t_avg,x,'r');

axes(ax(2));
plot(t_avg,y,'r');
