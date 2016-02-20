function spd = getLinSpd(cfg_in,pos)
% function spd = getLinSpd(cfg_in,pos)
%
% INPUT:
% pos tsd with x and y fields
%
% OUTPUT:
% tsd with spd field (linear speed in pixels/sec)
%
% cfg options with defaults:
%
% cfg.verbose = 1; 1 display command window text, 0 don't
%
% MvdM 2014-07-29

cfg_def.verbose = 1;
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose; fprintf('%s: calculating linear speed from position data\n',mfun); end

vx = dxdt(pos.tvec,getd(pos,'x'),'verbose',cfg.verbose);
vy = dxdt(pos.tvec,getd(pos,'y'),'verbose',cfg.verbose);

spd = tsd;
spd.tvec = pos.tvec;
spd.data = sqrt(vx.^2+vy.^2);
spd.label = {'spd'};

spd.cfg.history.mfun = cat(1,spd.cfg.history.mfun,mfun);
spd.cfg.history.cfg = cat(1,spd.cfg.history.cfg,{cfg});