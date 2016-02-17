function iv_out = DifferenceIV(cfg_in,iv1,iv2)
% function iv_out = DifferenceIV(cfg,iv1,iv2)
%
% keep only those iv1's that DO NOT include a piece of iv2
%
%  iv1    ____________    _____  ______     ____________   ______
%  iv2        ________    _____                _______             ______
%
%  iv_out                        ______                    ______
%
% INPUTS:
%
% iv1: interval data to be selected from based on..
% iv2: interval data to be used to select from iv1
%
% CFG OPTIONS:
% cfg.verbose = 1; If 1, displays text in command window; if 0, doesn't.
%
% OUTPUTS:
%
% iv_out: output interval data
%
% MvdM 2014-06-28

cfg_def.verbose = 1;
mfun = mfilename;

cfg = ProcessConfig(cfg_def,cfg_in,mfun); % should take whatever is in cfg_in and put it into cfg!

if ~CheckIV(iv1) || ~CheckIV(iv2); error('Inputs must be iv data type'); end

keep = ones(length(iv1.tstart),1);
for iI = 1:length(iv1.tstart)
    
    % first, check if any start or end is within this interval
    temp1 = find(iv2.tstart >= iv1.tstart(iI) & iv2.tstart <= iv1.tend(iI));
    temp2 = find(iv2.tend >= iv1.tstart(iI) & iv2.tend <= iv1.tend(iI));
    
    if ~isempty(temp1) | ~isempty(temp2)
       keep(iI) = 0;
       continue;
    end
    
    % check if interval is enveloped by anything
    for iJ = 1:length(iv2.tstart)
        if iv2.tstart(iJ) < iv1.tstart(iI) & iv2.tend(iJ) > iv1.tend(iI)
            keep(iI) = 0;
            break;
        end
    end
    
end

iv_out = iv1;
iv_out.tstart = iv_out.tstart(logical(keep));
iv_out.tend = iv_out.tend(logical(keep));

if cfg.verbose; fprintf('%s: %d intervals in, %d intervals out\n',mfun,length(iv1.tstart) + length(iv2.tstart),length(iv_out.tstart)); end

% housekeeping
iv1.cfg.history.mfun = cat(1,iv1.cfg.history.mfun,mfun);
iv1.cfg.history.cfg = cat(1,iv1.cfg.history.cfg,{cfg});