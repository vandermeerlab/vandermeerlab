function iv = addIV(cfg_in,iv)
% function iv = addIV(cfg_in,iv)
%
% expand/contract intervals by specified amount
%
% cfg.d: amount to expand/contract; must be two-element array which will be added to start and end iv fields respectively
%
% MvdM 2014-08-28 initial version


cfg_def = [];
cfg = ProcessConfig2(cfg_def,cfg_in);

mfun = mfilename;

if length(cfg.d) == 2
    iv.tstart = iv.tstart + cfg.d(1);
    iv.tend = iv.tend + cfg.d(2);
else
    error('cfg.d must contain 2 elements.');
end

iv_diffs = iv.tend-iv.tstart;
keep = find(iv_diffs > 0);

if length(keep) < length(iv.tend)
    warning(sprintf('%d/%d intervals discarded',length(iv.tend)-length(keep),length(iv.tend)));
    iv.tstart = iv.tstart(keep);
    iv.tend = iv.tend(keep);
end

% housekeeping
iv.cfg.history.mfun = cat(1,iv.cfg.history.mfun,mfun);
iv.cfg.history.cfg = cat(1,iv.cfg.history.cfg,{cfg});