function [ z, z_iv ] = linearize_shortcut(rat_id, unique_folder, pos_tsd, expkeys)
% Requires:
% - LinearizePos.m *Returns linearized position
% - makecoord_shortcut.m *Returns coord
% - restrict.m *Returns tsd data restricted to specific parts
%
% * Example usage:
% work_data = 'E:\data-shortcut\data-working\Shortcut-20150727';
% work_code = 'E:\code\shortcut';
% path_code = work_code;
% path_data = work_data;
% rat_id = 'R068_EI';
% cd(fullfile(path_code,'expkeys'))
% expday = expday(rat_id);
% unique_folder = expday.one;
% unique_id = unique_folder(1:15);
% expkeys = loadExpKeys_shortcut(unique_folder);
% cd(fullfile(path_data, rat_id, unique_folder));
% if exist([unique_id,'-vt.mat'],'file');
%     fprintf('*-vt.mat file found, loading.\n');
%     load([unique_id,'-vt.mat']);
% else
%     pos_tsd = position_shortcut(unique_folder,expkeys);
% end
% [z, z_iv] = linearize_shortcut(rat_id, unique_folder, pos_tsd, expkeys);
% 
% * Returns z.shortcut, z.u, z.novel, z_iv.shortcut, z_iv.u, z_iv.novel 
%   for shortcut experiment.

% Get conversion factor of pixel-to-cm from expkeys
conv_cm = expkeys.pxl_to_cm;

% Convert pos_tsd to cm
pos_tsd.data(1,:) = pos_tsd.data(1,:) ./ conv_cm(1);
pos_tsd.data(2,:) = pos_tsd.data(2,:) ./ conv_cm(2);

% Get trajectories from expkeys for each path
trajectory_pts(:,1) = expkeys.trajectory_pts(:,1) ./ conv_cm(1);
trajectory_pts(:,2) = expkeys.trajectory_pts(:,2) ./ conv_cm(2);
trajectory_labels = expkeys.trajectory_labels;
trajectory_idx = expkeys.trajectory_idx;

u_path = short_path(trajectory_pts,trajectory_labels,...
    trajectory_idx,'feeder1','feeder2');
shortcut_path = short_path(trajectory_pts,trajectory_labels,...
    trajectory_idx,'shortcut1','shortcut2');
novel_path = short_path(trajectory_pts,trajectory_labels,...
    trajectory_idx,'novel1','novel2');

% Set coodrinates for each trajectory (u, shortcut, novel)
coord_u = makecoord_shortcut(getd(pos_tsd,'x'),getd(pos_tsd,'y'),u_path);
coord_shortcut = makecoord_shortcut(getd(pos_tsd,'x'),getd(pos_tsd,'y'),...
    shortcut_path);
coord_novel = makecoord_shortcut(getd(pos_tsd,'x'),getd(pos_tsd,'y'),novel_path);

% Restrict to experiment Phase 3
phase3 = expkeys.phase3;
pos_tsd = restrict(pos_tsd, phase3(1), phase3(2));
pos_tsd = rm_nan(pos_tsd);

% Set distance threshold from "ideal" track
dist_thres = 7;

% Get linearized position and idx of points below distance threshold
cfg.Coord = coord_shortcut;
shortcut_lin = LinearizePos(cfg, pos_tsd); %get z_dist

cfg.Coord = coord_u;
u_lin = LinearizePos(cfg, pos_tsd);

cfg.Coord = coord_novel;
novel_lin = LinearizePos(cfg, pos_tsd);

cfg = [];
cfg.method = 'raw';
cfg.threshold = dist_thres;
cfg.dcn = '<';
cfg.target = 'z_dist';
z_iv.shortcut = TSDtoIV(cfg, shortcut_lin);
z_iv.u = TSDtoIV(cfg, u_lin);
z_iv.novel = TSDtoIV(cfg, novel_lin);

z.shortcut = restrict(shortcut_lin, z_iv.shortcut);
z.u = restrict(u_lin, z_iv.u);
z.novel = restrict(novel_lin, z_iv.novel);

save([unique_folder(1:15),'-linearized-xy.mat'], 'z', 'z_iv');
end

function [ all_path ] = short_path(trajectory_pts, trajectory_labels, ...
    trajectory_idx, start_label, end_label)
start_path = trajectory_idx(strcmp(trajectory_labels,start_label));
end_path = trajectory_idx(strcmp(trajectory_labels,end_label));
all_path = trajectory_pts(start_path:end_path, :);
end

function [ pos_tsd ] = rm_nan(pos_tsd)
keep = ~isnan(pos_tsd.data(1,:)) & ~isnan(pos_tsd.data(2,:));
pos_tsd.data = pos_tsd.data(:,keep);
pos_tsd.tvec = pos_tsd.tvec(keep);
end