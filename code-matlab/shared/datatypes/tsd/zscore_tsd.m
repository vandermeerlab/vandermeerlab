function tsd_in = zscore_tsd(tsd_in)
% function tsd_in = zscore_tsd(tsd_in)
%
%
% DOES NOT YET WORK FOR MULTIDIMENSIONAL DATA

cfg = [];
mfun = mfilename;

keep_idx = ~isnan(tsd_in.data);

tsd_in.data(keep_idx) = zscore(tsd_in.data(keep_idx));

% housekeeping
tsd_in.cfg.history.mfun = cat(1,tsd_in.cfg.history.mfun,mfun);
tsd_in.cfg.history.cfg = cat(1,tsd_in.cfg.history.cfg,{cfg});