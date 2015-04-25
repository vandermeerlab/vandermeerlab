% Requires:
% - event file (here '*Events.nev')
% - xy_targets_chunks(unique_folder,chunksize) (.m)
% - light_on(position_x,position_y,Timestamps,event_file,event_id) (.m)
% - leds_position(position_data,led1_position,led2_position,threshold) (.m)
% - make_tsd(pos_x,pos_y,Timestamps) (.m) *Not tested*

% *- position_constraint(position_data,window_size,outlier) (.m) *doesn't seem to help*     

clear all;

unique_folder = 'R068-2014-12-01_recording';
cd(['H:\data-working\Shortcut\R068_EI\',unique_folder]);

unique_id = unique_folder(1:15);
% position_file = sprintf('position-xy-%s',unique_id);
% load(position_file);

[pos_x,pos_y,Timestamps] = xy_targets_chunks(unique_folder,5000,400);

% Finding indices when lights are on (based on Event file)
light1_idx = light_on(pos_x,pos_y,Timestamps,'*Events.nev',10);
light2_idx = light_on(pos_x,pos_y,Timestamps,'*Events.nev',11);

% Converting lights from indices to a logical array
lights = false(size(pos_x));
lights(light1_idx) = true;
lights(light2_idx) = true;

% Finding LED constant position (*values may need to be altered based on exp*)
leds = leds_position(pos_x,pos_y,515,460,615,60,30);

% Finding outliers from position data
% outlier_x_idx = position_constraint(pos_x,4,5);
% outlier_y_idx = position_constraint(pos_y,4,5);

% Converting outliers from indices to a logical array
% outlier_x = false(size(pos_x));
% outlier_x(outlier_x_idx) = true;
% outlier_y = false(size(pos_y));
% outlier_y(outlier_y_idx) = true;

% Determining indices where the lights are on and in the LEDs locations.
% Removing annoying light points from analysis
annoying_x = leds & lights;
annoying_y = leds & lights;
annoying = annoying_x | annoying_y;

pos_x(annoying) = nan;
pos_y(annoying) = nan;

total_annoying = sum(annoying);
percent_annoying = ceil((sum(isnan(pos_x)) / length(pos_x)) * 100);


% Plotting to check
fig = figure('Position',[100, 100, 950, 950]);

subaxis(3,6,[2:5], 'Spacing', 0.04, 'Padding', 0.01, 'Margin', 0.04);
axis tight;
plot(pos_x,pos_y,'b.','MarkerSize',4);
pos_title = sprintf('Maze position of %s', unique_id);
set(gca,'xtick',[],'ytick',[]);
title(pos_title,'FontSize',14);

subaxis(3,1,2, 'Spacing', 0.04, 'Padding', 0.01, 'Margin', 0.04);
axis tight;
plot(Timestamps,pos_x,'k.','MarkerSize',4);
xpos_time_title = sprintf('X position of %s over time', unique_id);
title(xpos_time_title,'FontSize',14);

subaxis(3,1,3, 'Spacing', 0.04, 'Padding', 0.01, 'Margin', 0.04);
axis tight;
plot(Timestamps,pos_y,'k.','MarkerSize',4);
ypos_time_title = sprintf('Y position of %s over time', unique_id);
title(ypos_time_title,'FontSize',14);

% Make tsd
[pos_tsd,tvec] = make_tsd(pos_x,pos_y,Timestamps);


