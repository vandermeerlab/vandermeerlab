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
% (none yet)
%
% MvdM 2014-07-29

cfg_def = [];
cfg = ProcessConfig2(cfg_def,cfg_in);
mfun = mfilename;

vx = dxdt(pos.tvec,getd(pos,'x'));
vy = dxdt(pos.tvec,getd(pos,'y'));

spd = tsd;
spd.tvec = pos.tvec;
spd.data = sqrt(vx.^2+vy.^2);
spd.label = {'spd'};

spd.cfg.history.mfun = cat(1,spd.cfg.history.mfun,mfun);
spd.cfg.history.cfg = cat(1,spd.cfg.history.cfg,{cfg});