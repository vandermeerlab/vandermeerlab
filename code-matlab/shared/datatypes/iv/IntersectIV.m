function iv1 = IntersectIV(cfg_in,iv1,iv2)
% function iv = IntersectIV(cfg,iv1,iv2)
%
% keep only those iv1's that include a piece of iv2
%
% INPUTS:
%
% iv1: interval data to be selected from based on..
% iv2: interval data to be used to select from iv1
%
% CFG OPTIONS:
%
% OUTPUTS:
%
% iv_out: output interval data
%
% MvdM 2014-06-24

cfg_def = [];

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % should take whatever is in cfg_in and put it into cfg!

keep = zeros(length(iv1.tstart),1);
for iI = 1:length(iv1.tstart)
    
    % first, check if any start or end is within this interval
    temp1 = find(iv2.tstart > iv1.tstart(iI) & iv2.tstart < iv1.tend(iI));
    temp2 = find(iv2.tend > iv1.tstart(iI) & iv2.tend < iv1.tend(iI));
    
    if ~isempty(temp1) || ~isempty(temp2)
       keep(iI) = 1;
       continue;
    end
    
    % check if interval is enveloped by anything
    for iJ = 1:length(iv2.tstart)
        if iv2.tstart(iJ) < iv1.tstart(iI) && iv2.tend(iJ) > iv1.tend(iI)
            keep(iI) = 1;
            break;
        end
    end
    
end

iv1.tstart = iv1.tstart(logical(keep));
iv1.tend = iv1.tend(logical(keep));

% housekeeping
iv1.cfg.history.mfun = cat(1,iv1.cfg.history.mfun,mfun);
iv1.cfg.history.cfg = cat(1,iv1.cfg.history.cfg,{cfg});