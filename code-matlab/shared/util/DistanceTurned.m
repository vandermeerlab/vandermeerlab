function [distance_turned, prev_fd] = DistanceTurned(cfg_in, fd, fd_extra)
% [distance_turned, prev_fd] = DistanceTurned(cfg, fd, fd_extra)
%
% returns distance turned (in um) compared to previous session if available, Inf if otherwise
%
% required cfg fields:
%  cfg.this_rat_ID: numeric ID of current rat (as returned by getDataPath.m)
%  cfg.this_date: numeric ID of current session date (as returned by getDataPath.m)
%  cfg.ttno: tetrode number

cfg_def.this_rat_ID = [];
cfg_def.this_date = [];
cfg_def.ttno = [];
cfg_def.verbose = 0;

cfg = ProcessConfig(cfg_def, cfg_in);

LoadExpKeys;
today_depth = ExpKeys.TetrodeDepths(cfg.ttno);

idx_pool = 1:length(fd); % available fd idxs
idx_pool = idx_pool(ismember(idx_pool, find(fd_extra.ratID_num == cfg.this_rat_ID))); % filter by which rat we're looking at
idx_pool = idx_pool(ismember(idx_pool, find(fd_extra.fd_date_num < cfg.this_date))); % keep only preceding days
[~, prev_day_idx] = max(fd_extra.fd_date_num(idx_pool)); % find idx of closest day in pool

if isempty(prev_day_idx)
    if cfg.verbose
        fprintf('DistanceTurned: No previous session found.\n');
    end
    distance_turned = Inf; prev_fd = [];
    return;
end

prev_day_idx = idx_pool(prev_day_idx); % convert back to fd idx

if cfg.verbose
    fprintf('Previous session is %s\n',fd{prev_day_idx});
end

prev_fd = fd{prev_day_idx};
pushdir(prev_fd);
LoadExpKeys;
prev_depth = ExpKeys.TetrodeDepths(cfg.ttno);
popdir;

distance_turned = today_depth - prev_depth;