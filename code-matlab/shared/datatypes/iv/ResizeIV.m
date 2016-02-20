function iv = ResizeIV(cfg_in,iv)
% ResizeIV expand or contract intervals by specified amount
% function iv = ResizeIV(cfg_in,iv)
%
%   INPUTS
%      iv - iv struct containing interval data
% 
%   CONFIG options/defaults:
%
%      cfg.d = [0 0]; amount to expand/contract; must be two-element array 
%               which will be added to tstart and tend iv fields
%               respectively.
%               ex: [-1 1]  expand
%                   [1 -1]  contract
%
%               note that intervals can be shifted in cases like these:
%                   [-1 -1] shift left
%                   [1 1]   shift right 
% 
%      cfg.allowOverlap = 1; if 1, IV output may overlap if cfg.d is large; 
%               if 0, IVs will touch but not overlap
%
%         if cfg.allowOverlap is set to 1 the intervals might overlap:
%         tstart(1) ______________________tend(1)
%                       tstart(2)______________________tend(2)
%
%         if cfg.allowOverlap is set to 0 the intervals will touch:
%
%         tstart(1) __________________tend(1)
%                            tstart(2)_________________tend(2)
%
%
% MvdM 2014-08-28 initial version (originally called addIV)
% aacarey edit Mar 2015, added allowOverlap and merge
% aacarey edit Sept 2015, renamed addIV to ResizeIV and removed merge option

% make sure the input was created using the iv constructor
is_iv = CheckIV(iv);

if is_iv % proceed
    
    cfg_def.d = [0 0];
    cfg_def.allowOverlap = 1;
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
    
    % check/correct overlap
    if ~cfg.allowOverlap
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
end