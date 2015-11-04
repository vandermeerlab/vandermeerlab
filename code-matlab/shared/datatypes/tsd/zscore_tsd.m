function tsd_in = zscore_tsd(tsd_in)
% function tsd_out = zscore_tsd(tsd_in)
%
% z-scores data
%
% MvdM 2015-11-03 added functionality for multiple data channels

cfg = [];
mfun = mfilename;

nDim = size(tsd_in.data,1);

for iDim = 1:nDim
    
    keep_idx = ~isnan(tsd_in.data(iDim,:));
    
    tsd_in.data(iDim,keep_idx) = zscore(tsd_in.data(iDim,keep_idx));

end

% housekeeping
tsd_in.cfg.history.mfun = cat(1,tsd_in.cfg.history.mfun,mfun);
tsd_in.cfg.history.cfg = cat(1,tsd_in.cfg.history.cfg,{cfg});