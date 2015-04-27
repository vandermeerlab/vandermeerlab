%%
cd('D:\data\R066\R066-2014-11-27');

%% load raw data
cfg = [];
cfg.targets = 1:4;
cfg.mode = 'raw';

[ x_avg,y_avg,t_avg ] = xy_targets_chunks(cfg);
t_avg = t_avg * 10^-6; % convert to s

%% plot each target individually (raw)
cols = 'rgbcm';
for iNT = 4:-1:1
    plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT));
    hold on;
end

%% find intervals with LEDs on

cfg = []; 
cfg.eventList = {'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0001).', ...
    'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0002).', ...
    'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0000).'};

evt = LoadEvents(cfg);

led0 = iv(evt.t{1},evt.t{3}(nearest_idx3(evt.t{1},evt.t{3},1)));
led1 = iv(evt.t{2},evt.t{3}(nearest_idx3(evt.t{2},evt.t{3},1)));

%% check
subplot(121)
hist(led0.tend-led0.tstart,100)
subplot(122)
hist(led1.tend-led1.tstart,100)

%% expand IVs -- this is apparently needed to catch stray position samples which fall outside the LED event window
cfg = []; cfg.d = [-1/90 1/90];
led0 = addIV(cfg,led0);
led1 = addIV(cfg,led1);

%% restrict
xtsd = tsd(t_avg,x_avg);
led0_idx = TSD_getidx(xtsd,led0);
led1_idx = TSD_getidx(xtsd,led1);

%% plot
plot(t_avg(led0_idx),x_avg(led0_idx,1),'o','Color',cols(1));
plot(t_avg(led1_idx),x_avg(led1_idx,1),'o','Color',cols(1));

%% define LED centroids (take from plots)
led0_xy_centroid = [];
led1_xy_centroid = [];

%% get distance of each sample to LEDs

%% if distance below threshold (radius) for either LED, then keep

%% set position sample to nan if LED is on AND position is within radius
