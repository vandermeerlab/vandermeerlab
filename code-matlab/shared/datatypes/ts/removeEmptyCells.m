function [S,keep_idx] = removeEmptyCells(S)
% function [S,keep_idx] = removeEmptyCells(S)

cfg = [];
mfun = mfilename;

keep_idx = ~cellfun(@isempty,S.t);
S.t = S.t(keep_idx);
S.label = S.label(keep_idx);

% housekeeping
S.cfg.history.mfun = cat(1,S.cfg.history.mfun,mfun);
S.cfg.history.cfg = cat(1,S.cfg.history.cfg,{cfg});
