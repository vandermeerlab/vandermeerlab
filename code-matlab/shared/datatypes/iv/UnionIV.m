function iv1 = UnionIV(cfg_in,iv1,iv2)
% function iv = UnionIV(cfg,iv1,iv2)
%
% union of iv objects
%
% output is resorted based on start times (ascending)
%
% MvdM 2014-08-28 initial version

cfg_def = [];
cfg = ProcessConfig(cfg_def,cfg_in);

mfun = mfilename;

if isempty(iv1) % function should work for empty arguments
    iv1 = iv([],[]);
end

if isempty(iv2)
    iv2 = iv([],[]);
end

if ~CheckIV(iv1) | ~CheckIV(iv2)
   error('Malformed IV.'); 
end

iv1.tstart = cat(1,iv1.tstart,iv2.tstart); % note that constructor guarantees column vectors
iv1.tend = cat(1,iv1.tend,iv2.tend);

[iv1.tstart,sort_idx] = sort(iv1.tstart,'ascend');
iv1.tend = iv1.tend(sort_idx);

% housekeeping
iv1.cfg.history.mfun = cat(1,iv1.cfg.history.mfun,mfun);
iv1.cfg.history.cfg = cat(1,iv1.cfg.history.cfg,{cfg});