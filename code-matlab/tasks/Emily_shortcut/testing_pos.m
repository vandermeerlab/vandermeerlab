clear all;
rat_id = 'R063_EI';
expday = expday(rat_id);
unique_folder = expday.eight;
unique_id = unique_folder(1:15);
ExpKeys = load_ExpKeys(rat_id,unique_folder);

cd(fullfile('H:\data-working\Shortcut',rat_id,unique_folder));

cfg = [];
cfg.write_output = 1;
cfg.led0_xy = ExpKeys.led1_xy;
cfg.led1_xy = ExpKeys.led2_xy;
cfg.thr = ExpKeys.led_radius;
vt = ExtractShortcutVT(cfg);

% Plotting to check
fig = figure();
width = 10;
height = 10;
fontsize = 15;
linewidth = 1.5;
markersize = 4;

set(0,'defaultLineLineWidth',linewidth);
set(0,'defaultLineMarkerSize',markersize);
set(0,'defaultLineLineWidth',linewidth);
set(0,'defaultLineMarkerSize',markersize);

defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [100, 100, 950, 950]);

set(0,'defaultFigureInvertHardcopy','on');
set(0,'defaultFigurePaperUnits','inches');
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

subaxis(3,8,[3:6], 'Spacing', 0.04, 'Padding', 0.04, 'Margin', 0.04);
axis tight;
plot(vt.data(1,:),vt.data(2,:),'b.','MarkerSize',markersize);
pos_title = sprintf('Maze position of %s', unique_id);
set(gca,'xtick',[],'ytick',[]);
title(pos_title,'FontSize',fontsize);

subaxis(3,1,2, 'Spacing', 0.04, 'Padding', 0.04, 'Margin', 0.04);
axis tight;
plot(vt.tvec,vt.data(1,:),'k.','MarkerSize',markersize);
xpos_time_title = sprintf('X position of %s over time', unique_id);
title(xpos_time_title,'FontSize',fontsize);
xlabel('time (s)');
ylabel('position (au)');

subaxis(3,1,3, 'Spacing', 0.04, 'Padding', 0.04, 'Margin', 0.04);
axis tight;
plot(vt.tvec,vt.data(2,:),'k.','MarkerSize',markersize);
ypos_time_title = sprintf('Y position of %s over time', unique_id);
title(ypos_time_title,'FontSize',fontsize);
xlabel('time (s)');
ylabel('position (au)');
