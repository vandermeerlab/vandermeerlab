function [ pos_tsd ] = emi_position(unique_folder,expkeys)
% Requires:
% - emi_lights.m *Returns event_times
%
% * Example usage:
% rat_id = 'R068_EI';
% expday = emi_expday(rat_id);
% unique_folder = expday.one;
% unique_id = unique_folder(1:15);
% cd(fullfile('C:\Users\Emily\Desktop',rat_id,unique_folder));
% ExpKeys = loadExpKeys_shortcut(unique_folder);
% pos_tsd = position_shortcut(unique_folder,ExpKeys);
% * Returns and saves pos_tsd for shortcut experiment
unique_id = unique_folder(1:15);

num_target = 2;
led1_xy = expkeys.led1_xy;
led2_xy = expkeys.led2_xy;
led_radius = expkeys.led_radius;
head_dist = 50;
light_interval = 0.1; % from nlx initializemaze
cushion = 0.01; % for evt times & evt intervals

if exist(['pos-xytime-',unique_id,'.mat'],'file');
    fprintf('pos-xytime-*.mat file found, loading. \n');
    load(['pos-xytime-',unique_id,'.mat']);
else
    [pos_x, pos_y, Timestamps] = xy_targets_shortcut(unique_folder,num_target);
end

[led1on_times, led2on_times, ledoff_times] = light_on(unique_folder, expkeys);

led1off_times = ledoff_times(nearest_idx3(led1on_times, ledoff_times));
led1iv = iv(led1on_times, led1off_times);
led1iv_corrected = false(1, length(led1iv.tstart));
led1iv_corrected((led1iv.tend - led1iv.tstart) >=  light_interval-cushion) = true;
led1iv = iv(led1on_times(led1iv_corrected)-cushion, led1off_times(led1iv_corrected)+cushion);

led2off_times = ledoff_times(nearest_idx3(led2on_times,ledoff_times));
led2iv = iv(led2on_times, led2off_times);
led2iv_corrected = find((led2iv.tend - led2iv.tstart) >=  light_interval-cushion);
led2iv = iv(led2on_times(led2iv_corrected)-cushion, led2off_times(led2iv_corrected)+cushion);

% Get idxs corresponding to "LED on"
xtsd = tsd(Timestamps, pos_x);

led1_keep = false(size(xtsd.tvec));
for u1 = 1:length(led1iv.tstart)
    led1_keep = led1_keep | (xtsd.tvec >= led1iv.tstart(u1) & xtsd.tvec <= led1iv.tend(u1));
end
lights1 = led1_keep;

led2_keep = false(size(xtsd.tvec));
for u2 = 1:length(led2iv.tstart)
    led2_keep = led2_keep | (xtsd.tvec >= led2iv.tstart(u2) & xtsd.tvec <= led2iv.tend(u2));
end
lights2 = led2_keep;

% Get euclidian distance of each position sample to LED reward center
DistFun_posled = @(xd,yd,ctr) sqrt((xd-ctr(1)).^2 + (yd-ctr(2)).^2);

led1_dist = nan(num_target, length(pos_x));
led2_dist = nan(num_target, length(pos_x));
for i = 1:num_target
    led1_dist(i,:) = DistFun_posled(pos_x(i,:), pos_y(i,:), led1_xy);
    led2_dist(i,:) = DistFun_posled(pos_x(i,:), pos_y(i,:), led2_xy);
end

% Logical array for samples within LED position
leds1 = false(size(pos_x));
leds1(led1_dist <= led_radius) = true;
leds2 = false(size(pos_x));
leds2(led2_dist <= led_radius) = true;

% Set position data to nan is LED is on & position is within led location
for j = 1:num_target
    nan_idx = (lights1 & leds1(j,:)) | (lights2 & leds2(j,:));
    pos_x(j,nan_idx) = nan;
    pos_y(j,nan_idx) = nan;
end

% Set position date to nan when jump & position is within led location
finish_removing = false;

while finish_removing == false
    removed = 0;
    for k = 1:num_target
        continuity = cat(2,0,sqrt(diffskipnan(pos_x(k,:)).^2 + ...
            diffskipnan(pos_y(k,:)).^2));
        continuity_idx = (continuity >= head_dist) & (leds1(k,:)|leds2(k,:));
        pos_x(k,continuity_idx) = nan;
        pos_y(k,continuity_idx) = nan;
        
        removed = removed + sum(continuity_idx);
    end
    
    if (removed == 0)
        finish_removing = true;
    end
end

% Average between targets when close; otherwise take 1st target
avg_x = nan(1,length(pos_x));
avg_y = nan(1,length(pos_y));

if num_target == 1 % Not tested...
    avg_x = pos_x;
    avg_y = pos_y;
elseif num_target == 2
    DistFun_targ = @(targ1,targ2) sqrt((targ1(1,:)-targ2(1,:)).^2 + ...
        (targ1(2,:)-targ2(2,:)).^2);
    distance = DistFun_targ([pos_x(1,:); pos_y(1,:)], [pos_x(2,:); pos_y(2,:)]);
    
    far_idx = distance >= head_dist;
    avg_x(far_idx) = pos_x(1,far_idx);
    avg_y(far_idx) = pos_y(1,far_idx);

    other_idx = ~far_idx;
    avg_x(other_idx) = nanmean(pos_x(:,other_idx));
    avg_y(other_idx) = nanmean(pos_y(:,other_idx));
else
    error('Needs 1 or 2 targets. Only one target is not tested.');
end

% Filtering with median filter to remove stragglers (and put back nans)
only_nan = isnan(avg_x);
filter_radius = 6;
avg_x = medfilt1m(avg_x, filter_radius);
avg_y = medfilt1m(avg_y, filter_radius);
avg_x(only_nan) = nan; 
avg_y(only_nan) = nan;

pos_tsd = pos_tsd_shortcut(avg_x, avg_y, Timestamps);

if ~CheckTSD(pos_tsd)
    error('tsd check failed.');
end

save([unique_folder(1:15),'-vt.mat'],'pos_tsd');
end

function [ pos_tsd ] = pos_tsd_shortcut( pos_x,pos_y,Timestamps )
% Uses datatype tsd.m function
pos_tsd = tsd;
pos_tsd.tvec = Timestamps;
pos_tsd.data(1,:) = pos_x;
pos_tsd.label{1} = 'x';
pos_tsd.data(2,:) = pos_y;
pos_tsd.label{2} = 'y';
end
