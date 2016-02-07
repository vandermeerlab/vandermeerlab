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

% Get positions
if exist([unique_id,'-vt.mat'],'file');
    fprintf('*-vt.mat file found, loading.\n');
    load([unique_id,'-vt.mat']);
else
    pos_tsd = position_shortcut(unique_folder,expkeys);
end

% Restrict data to Phase 3
pos_tsd = restrict(pos_tsd, expkeys.phase3(1), expkeys.phase3(2));

% Get linearized "ideal" position
[z, z_iv] = linearize_shortcut(rat_id, unique_folder, pos_tsd, expkeys);

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

% Get ordered spikes for plotting
spikes_sorted = spikes_filtered;
spikes_sorted.t = spikes_sorted.t(tc_sort_idx);
spikes_sorted.label = spikes_sorted.label(tc_sort_idx);

% Quantitative raster with MakeQfromS.m
dt = median(diff(zlin.tvec));
cfg = [];
cfg.tvec_edges = [zlin.tvec-dt./2,zlin.tvec(end)];
subsample = 6;
cfg.tvec_edges = cfg.tvec_edges(1:subsample:end);
cfg.smooth = 'gauss';
raster = MakeQfromS(cfg, spikes_filtered);

% Bayesian decoder, using DecodeZ()
cfg = [];
cfg.noSpikesInBin = 'nans';
decoded_norm = DecodeZ(cfg, raster, tc.tc);

[decoded, max_decoded_idx] = max(decoded_norm.data, [], 2);

[z_binned, z_binned_idx] = histc(zlin.data, z_min:binsize:z_max, 6);        
[z_binned, z_binned_idx] = trim_histc(z_binned, z_binned_idx);

firing = ~isnan(decoded);
decoded_idx = max_decoded_idx(firing);
z_idx = z_binned_idx(firing);

decode_error = abs(decoded_idx'-z_idx);
decode_error_avg = nanmean(decode_error);

decoded_norm.data(isnan(decoded_norm.data)) = 0;


figure(1); clf; hold on;
py_decoded = decoded_norm.data';
imagesc(decoded_norm.data');
xlim([0, size(decoded_norm.data,1)]);
ylim([0, size(decoded_norm.data,2)]);
set(gca,'FontSize',14); xlabel('Time (s)'); ylabel('Linearized position');
colormap(jet);


figure(2); clf; hold on;
imagesc(phase_start:dt:phase_end,...
 1:size(raster.data(tc_sort_idx,:),1),raster.data(tc_sort_idx,:));
py_raster = raster.data(tc_sort_idx,:);
colormap(jet);
ylim([0, size(raster.data(tc_sort_idx,:),1)]);
xlim([phase_start, phase_end]);
title('Quantitative raster of neurons along trajectory');


figure(3); clf; hold on;
plot(raster.tvec(firing), z_idx, 'k.','markersize',10);
plot(decoded_norm.tvec(firing), decoded_idx, 'r.','markersize',5);
title('Decoded position. Testing decoder.');


figure(4); clf; hold on; 
% ylim([-33, length(tc_max_idx(tc_sort_idx))]);
% ylabel([0, length(tc_max_idx(tc_sort_idx)), 20]);
pos_max = 15;
zlin.data = rescale(zlin.data, -pos_max, 0);
plot(zlin.tvec, zlin.data, 'r*','markersize', 2);
axis off;
decoded_idx = rescale(decoded_idx, -pos_max+2, -pos_max+17);
plot(decoded_norm.tvec(firing), decoded_idx, 'b*', 'markersize', 2);
xlim([26930, 27200]);
ylim([-17, 3]);
% legend([zlin.data,decoded_idx],'Actual position','Decoded position');


