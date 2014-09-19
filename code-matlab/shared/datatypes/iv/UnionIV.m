function iv1 = UnionIV(cfg_in,iv1,iv2)
% function iv = UnionIV(cfg,iv1,iv2)
%
% union of iv objects
%
% output is resorted based on start times (ascending)
%
% MvdM 2014-08-28 initial version
% EC 2014-09-28 
%      Added a correction incase the dimesions were not correct
%      for the cat(1,...) command.  If this is the case then iv1 and iv2 
%      are transposed to allow for veritcal concatination.  

cfg = [];
ProcessConfig;

mfun = mfilename;
if length(iv1.tstart)==1 && length(iv2.tstart)==1
    iv1.tstart = cat(1,iv1.tstart',iv2.tstart');
    iv1.tend = cat(1,iv1.tend',iv2.tend');
else
    
    iv1.tstart = cat(1,iv1.tstart,iv2.tstart);
    iv1.tend = cat(1,iv1.tend,iv2.tend);
end
[iv1.tstart,sort_idx] = sort(iv1.tstart,'ascend');
iv1.tend = iv1.tend(sort_idx);

% housekeeping
iv1.cfg.history.mfun = cat(1,iv1.cfg.history.mfun,mfun);
iv1.cfg.history.cfg = cat(1,iv1.cfg.history.cfg,{cfg});