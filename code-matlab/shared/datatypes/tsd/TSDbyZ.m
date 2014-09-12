function z_av = TSDbyZ(cfg_in,z,tsd_in)
% function z_av = TSDbyZ(cfg,z,tsd_in)
%
% compute averages of the values in tsd_in binned by z's cfg.edges
%
% inputs:
%
% outputs:
%
% cfg options:
%
%

cfg = [];
cfg.edges = [];

ProcessConfig;

[~,z_idx] = histc(z.data(1,:),cfg.edges); % need to fix points that fall on last edge

tsd_z = interp1(tsd_in.tvec,tsd_in.data(1,:),z.tvec,'nearest'); % tsd value for each z sample

z_av = nan(length(cfg.edges)-1,1);

idx_list = unique(z_idx);
for iI = 1:length(idx_list)
    curr_idx = idx_list(iI);
    z_av(curr_idx) = nanmean(tsd_z(z_idx == curr_idx));
end


 