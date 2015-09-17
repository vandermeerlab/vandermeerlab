laptop_data = 'H:\data-working\Shortcut-20150611';
work_data = 'E:\data-shortcut\data-working\Shortcut-20150727';
laptop_code = 'C:\Users\Emily\Code\Shortcut';
work_code = 'E:\code\shortcut';
path_code = work_code;
path_data = work_data;

rat_id = 'R068_EI';
cd(fullfile(path_code,'expkeys'))
expday = emi_expday(rat_id);
unique_folder = expday.one;
unique_id = unique_folder(1:15);

expkeys = loadExpKeys_shortcut(unique_folder);

cd(fullfile(path_data, rat_id, unique_folder));

% Get position and restrict to phase3
if exist([unique_id,'-vt.mat'],'file');
    fprintf('*-vt.mat file found, loading.\n');
    load([unique_id,'-vt.mat']);
else
    pos_tsd = position_shortcut(unique_folder,expkeys);
end
pos_tsd = restrict(pos_tsd, expkeys.phase3(1), expkeys.phase3(2));

% Get linearized "ideal" position
if exist([unique_id,'-linearized-xy.mat'],'file');
    fprintf('*-linearized-xy.mat file found, loading.\n');
    load([unique_id,'-linearized-xy.mat']);
else
    [z, z_iv] = linearize_shortcut(rat_id, unique_folder, pos_tsd, expkeys);
end
% Tuning is confused by z_dist
z.u.data = z.u.data(1,:);
z.shortcut.data = z.shortcut.data(1,:);
z.novel.data = z.novel.data(1,:);

% Investigating U trials only for now
zlin = z.u;
zlin_iv = z_iv.u;

% Filtering neurons using spikes_filtered_shortcut.m
cfg = [];
cfg.load_questionable_cells = 1;
cfg.useClustersFile = 0;
spikes = LoadSpikes(cfg);
spikes_filtered = spikes_filtered_shortcut(spikes, pos_tsd, zlin, zlin_iv, expkeys);

% Tuning curves with TuningCurves.m
z_min = min(zlin.data(1,:));
z_max = max(zlin.data(1,:));
binsize = 3;
cfg = [];
cfg.binEdges{1} = z_min:binsize:z_max;
cfg.smoothingKernel = gausskernel(15,3); % (wide by std dev)
clear tc;
tc = TuningCurves(cfg, spikes_filtered, zlin);

% Sorting neurons based on peak firing along position
[~, tc_max_idx] = max(tc.tc, [], 2);
[~, tc_sort_idx] = sort(tc_max_idx);
sorted_tc = tc.tc(tc_sort_idx,:);

% Quantitative raster with MakeQfromS.m
dt = median(diff(zlin.tvec));
cfg = [];
cfg.tvec_edges = [zlin.tvec-dt./2,zlin.tvec(end)];
%cfg.tvec_edges = [cfg.tvec_edges(1:2:end),z.tvec(end)];
cfg.smooth = 'gauss';
raster = MakeQfromS(cfg, spikes_filtered);

% Plot to check
figure(1); clf; hold on;
imagesc(expkeys.phase3(1):dt:expkeys.phase3(end),...
    1:size(raster.data(tc_sort_idx,:),1),raster.data(tc_sort_idx,:));
colormap(jet);
ylim([0, size(raster.data(tc_sort_idx,:),1)]);
xlim([expkeys.phase3(1), expkeys.phase3(2)]);
title('Quantitative raster of neurons along U. Testing decoder.');

figure(2); clf; hold on;
plot(sorted_tc(4,:), 'b', 'LineWidth', 2);
plot(sorted_tc(34,:), 'g', 'LineWidth', 2);
plot(sorted_tc(59,:), 'm', 'LineWidth', 2);
title('Tuning curve of three neurons along U. Testing decoder.');


