function iv = addIV(cfg_in,iv)
% ADDIV expand/contract intervals by specified amount
% function iv = addIV(cfg_in,iv)
%
%   INPUTS
%      iv - iv struct containing interval data
% 
%   CONFIG options/defaults:
%
%      cfg.d = [0 0]; amount to expand/contract; must be two-element array 
%               which will be added to start and end iv fields respectively
%               ex: [-0.01 0.01] 
% 
%      cfg.allowOverlap = 1; if 1, IV output may overlap if cfg.d is large; 
%               if 0, IVs will touch but not overlap
%
%      cfg.merge = 0; if 0, events are not merged if they touch or overlap; 
%               if 1, events are merged (has precedence over allowOverlap)
%               Note: merging can be done in a separate step using
%               MergeSingleIV().
%
% MvdM 2014-08-28 initial version
% ACarey edit Mar 2015, added allowOverlap and merge


cfg_def.d = [0 0];
cfg_def.allowOverlap = 1;
cfg_def.merge = 0;
cfg = ProcessConfig2(cfg_def,cfg_in);

mfun = mfilename;

if length(cfg.d) == 2
    iv.tstart = iv.tstart + cfg.d(1);
    iv.tend = iv.tend + cfg.d(2);
else
    error('cfg.d must contain 2 elements.');
end

% remove "bad" intervals 
iv_diffs = iv.tend-iv.tstart;
keep = find(iv_diffs > 0);

if length(keep) < length(iv.tend)
    warning('%d/%d intervals discarded',length(iv.tend)-length(keep),length(iv.tend));
    iv.tstart = iv.tstart(keep);
    iv.tend = iv.tend(keep);
end

% merge if user wants
if cfg.merge
    iv = MergeSingleIV([],iv);
end

% check/correct overlap
if ~cfg.merge && ~cfg.allowOverlap
    for ii = 1:length(iv.tstart)-1
        if iv.tend(ii) > iv.tstart(ii+1)
            avg = mean([iv.tend(ii) iv.tstart(ii+1)]);
            iv.tend(ii) = avg;
            iv.tstart(ii+1) = avg;
        end
    end
end

% housekeeping
iv.cfg.history.mfun = cat(1,iv.cfg.history.mfun,mfun);
iv.cfg.history.cfg = cat(1,iv.cfg.history.cfg,{cfg});