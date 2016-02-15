function space_av = TSDbySpace(cfg_in,pos,tsd_in)
% function space_av = TSDbySpace(cfg,pos,tsd_in)
%
% construct 2D matrix with spatial averages of the values in tsd_in
%
% inputs:
%
% outputs:
%
% cfg options:
%
%

cfg_def = [];
cfg_def.x_edges = [];
cfg_def.y_edges = [];

cfg = ProcessConfig(cfg_def,cfg_in);

pos_mat(:,1) = getd(pos,'y');
pos_mat(:,2) = getd(pos,'x'); 

[occ_hist,edges,mid,pos_idx] = histcn(pos_mat,cfg.y_edges,cfg.x_edges);

tsd_pos = interp1(tsd_in.tvec,tsd_in.data,pos.tvec,'nearest'); % tsd value for each pos sample

space_av = nan(length(cfg.y_edges)-1,length(cfg.x_edges)-1); % need to fix points that fall on last edge though

temp = sub2ind(size(space_av),pos_idx(:,1),pos_idx(:,2)); % idx into binned space for each pos sample

idx_list = unique(temp);
for iI = 1:length(idx_list)
    curr_idx = idx_list(iI);
    space_av(curr_idx) = nanmean(tsd_pos(find(temp == curr_idx)));
end


 