function iv = SelectIV_idx(iv,keep_idx)
% function iv_out = SelectIV_idx(iv_inmkeep_idx)
%
% select IVs based idxs


iv.tstart = iv.tstart(keep_idx);
iv.tend = iv.tend(keep_idx);

if isfield(iv,'usr')
    for iU = 1:length(iv.usr)
        iv.usr(iU).data = iv.usr(iU).data(keep_idx);
    end
end
