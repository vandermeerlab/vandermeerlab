function [ pos ] = ExtractShortcutVT(cfg_in)
% function ExtractShortcutVT(cfg_in)
%
% extract VT data from Emily's shortcut task raw data
%
% run from data folder; assumes VT1.nvt and *Events.Nev present
%
% cfg_def.nTargets = 2; % number of videotracker targets to extract; 2 is usually the red and green rat LEDs
% cfg_def.plot = 0; % plot intermediate processing step results for debugging
% cfg_def.led0_xy = [521 473]; % this is ok for R066-2014-11-27
% cfg_def.led1_xy = [623 57];
% cfg_def.thr = 10; % distance (pixels) from LED coordinates to be considered a LED target
% cfg_def.jump_thr = 50; % distance (pixels) between position samples to be considered a LED target
% cfg_def.write_output = 1;
% cfg_def.dist = 20; % averaging threshold for distance between red and green targets
%
% MvdM 2015-05-03 initial version

cfg_def.nTargets = 2; % number of videotracker targets to extract; 2 is usually the red and green rat LEDs
cfg_def.plot = 0; % plot intermediate processing step results for debugging
cfg_def.led0_xy = [521 473]; % this is ok for R066-2014-11-27
cfg_def.led1_xy = [623 57];
cfg_def.thr = 10; % distance (pixels) from LED coordinates to be considered a LED target
cfg_def.jump_thr = 50; % distance (pixels) between position samples to be considered a LED target
cfg_def.write_output = 1;
cfg_def.dist = 20; % averaging threshold for distance between red and green targets

cfg = ProcessConfig2(cfg_def,cfg_in);


% 
if isempty(FindFile('*VT1.nvt'))
   error('No VT1.nvt file found.'); 
end

if isempty(FindFile('*events.nev'))
    error('No *events.nev file found.');
end

%
[~,fd,~] = fileparts(pwd);

%
cfg_loading = [];
cfg_loading.targets = 1:cfg.nTargets;
cfg_loading.mode = 'raw';

if exist('temp_pos.mat','file');
    fprintf('temp_pos.mat file found, loading...\n');
    load temp_pos;
else
    [x_avg,y_avg,t_avg] = xy_targets_chunks(cfg_loading);
    t_avg = t_avg * 10^-6; % convert to s
end

% plot each target individually (raw)
if cfg.plot
    
    figure(1); cols = 'rgbcm';
    
    subplot(221);
    for iNT = nTargets:-1:1
        plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT)); hold on;
    end; title('x (raw)');
    
    subplot(222);
    for iNT = nTargets:-1:1
        plot(t_avg,y_avg(:,iNT),'.','Color',cols(iNT)); hold on;
    end; title('y (raw)');
   
    subplot(223); hist(x_avg,1:720);
    subplot(224); hist(y_avg,1:640);
      
end

% find intervals with LEDs on
cfg_evt = []; 
cfg_evt.eventList = {'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0001).', ...
    'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0002).', ...
    'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0000).'};

evt = LoadEvents(cfg_evt);

led0 = iv(evt.t{1},evt.t{3}(nearest_idx3(evt.t{1},evt.t{3},1)));
led1 = iv(evt.t{2},evt.t{3}(nearest_idx3(evt.t{2},evt.t{3},1)));

if any(led0.tend-led0.tstart < 0) | any(led1.tend-led1.tstart < 0)
    warning('Negative LED intervals found.');
end

% expand IVs -- this is apparently needed to catch stray position samples which fall outside the LED event window
cfg_iv = []; cfg_iv.d = [-1/60 1/60];
led0 = ResizeIV(cfg_iv,led0); led1 = ResizeIV(cfg_iv,led1);

% get idxs corresponding to "LED on"
fprintf('Getting idxs for LED on times...\n');
xtsd = tsd(t_avg,x_avg);
led0_idx = TSD_getidx(xtsd,led0); led1_idx = TSD_getidx(xtsd,led1);

% plot each target individually (raw) plus points to be removed
if cfg.plot
    figure(2);
    ax(1) = subplot(211);
    for iNT = cfg.nTargets:-1:1
        plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT)); hold on;
    end
    plot(t_avg(led0_idx),x_avg(led0_idx,1),'o','Color',cols(1));
    plot(t_avg(led1_idx),x_avg(led1_idx,1),'o','Color',cols(1));
    
    ax(2) = subplot(212);
    for iNT = cfg.nTargets:-1:1
        plot(t_avg,y_avg(:,iNT),'.','Color',cols(iNT)); hold on;
    end
    plot(t_avg(led0_idx),y_avg(led0_idx,1),'o','Color',cols(1));
    plot(t_avg(led1_idx),y_avg(led1_idx,1),'o','Color',cols(1));
    
    linkaxes(ax,'x');
end


% get distance of each sample to LEDs
DistFun = @(xd,yd,ctr) sqrt((xd-ctr(1)).^2 + (yd-ctr(2)).^2); % euclidian distance is only accurate if pixels are square (which they are not)

for iT = cfg.nTargets:-1:1
    led0_dist(:,iT) = DistFun(x_avg(:,iT),y_avg(:,iT),cfg.led0_xy);
    led1_dist(:,iT) = DistFun(x_avg(:,iT),y_avg(:,iT),cfg.led1_xy);
end

% set position sample to nan if LED is on AND position is within radius
led0_bin = zeros(size(led0_dist(:,1))); led0_bin(led0_idx) = 1; % binarize for easier comparison
led1_bin = zeros(size(led1_dist(:,1))); led1_bin(led1_idx) = 1;

for iT = cfg.nTargets:-1:1
    
    nan_idx(:,iT) = (led0_bin & led0_dist(:,iT) < cfg.thr) | (led1_bin & led1_dist(:,iT) < cfg.thr);
    x_avg(nan_idx(:,iT),iT) = NaN; y_avg(nan_idx(:,iT),iT) = NaN; 

end

% remove samples which jump + are in LED locations
all_removed = 0; iI = 1;

while ~all_removed
    
    fprintf('Removing jumps, iteration %d...\n',iI);
    
    clear nan_idx dist_diff;
    for iNT = 2:-1:1
        % first get distances between successive samples (diffs)
        dist_diff(:,iNT) = cat(1,0,sqrt(diffskipnan(x_avg(:,iNT)).^2 + diffskipnan(y_avg(:,iNT)).^2));
        
        nan_idx(:,iNT) = (dist_diff(:,iNT) > cfg.jump_thr) & ((led0_dist(:,iNT) < cfg.thr) | (led1_dist(:,iNT) < cfg.thr));
        
        x_avg(nan_idx(:,iNT),iNT) = NaN; y_avg(nan_idx(:,iNT),iNT) = NaN;
        
    end
    
    if all(nan_idx(:) == 0)
        all_removed = 1;
    end
    iI = iI + 1;
end

% combine samples from red and green targets
%
% if two targets available and are plausibly close, then average target
% positions
x = nan(size(x_avg(:,1))); y = nan(size(y_avg(:,1))); % initialize output vars

DistFun2 = @(led1,led2) sqrt((led1(:,1)-led2(:,1)).^2 + (led1(:,2)-led2(:,2)).^2);
d = DistFun2(cat(2,x_avg(:,1),y_avg(:,1)),cat(2,x_avg(:,2),y_avg(:,2)));
avg_idx = d < cfg.dist;

x(avg_idx) = mean(cat(2,x_avg(avg_idx,1),x_avg(avg_idx,2)),2);
y(avg_idx) = mean(cat(2,y_avg(avg_idx,1),y_avg(avg_idx,2)),2);

% if targets not close, take red target
non_avg_idx = d > cfg.dist;
x(non_avg_idx) = x_avg(non_avg_idx,1); y(non_avg_idx) = y_avg(non_avg_idx,1);

% if one target available %(because of LED on) take that target
nnan = isnan(x_avg(:,1)) + isnan(x_avg(:,2));
led1t = (nnan == 1);% & (led0_bin | led1_bin);
x(led1t) = nanmean(cat(2,x_avg(led1t,1),x_avg(led1t,2)),2);
y(led1t) = nanmean(cat(2,y_avg(led1t,1),y_avg(led1t,2)),2);

% filter
nan_idx = isnan(x);

x = medfilt1m(x,6);
y = medfilt1m(y,6);

% put NaNs back
x(nan_idx) = nan; y(nan_idx) = nan;

if cfg.plot
    figure(3);
    
    ax(1) = subplot(211);
    for iNT = 2:-1:1
        plot(t_avg,x_avg(:,iNT),'.','Color',cols(iNT));
        hold on;
    end
    
    ax(2) = subplot(212);
    for iNT = 2:-1:1
        plot(t_avg,y_avg(:,iNT),'.','Color',cols(iNT));
        hold on;
    end
    
    axes(ax(1));
    plot(t_avg,x,'r');
    
    axes(ax(2));
    plot(t_avg,y,'r');
    
    linkaxes(ax,'x');
end

%
data(:,1) = x; data(:,2) = y; label = {'x','y'};
pos = tsd(t_avg,data',label);

if ~CheckTSD(pos)
   error('TSD check NOT passed!'); 
end

if cfg.write_output
    f_out = cat(2,fd,'-vt.mat');
    save(f_out,'pos');
end