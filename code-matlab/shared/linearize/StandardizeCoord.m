function Coord_out = StandardizeCoord(cfg_in,Coord_in)
% function Coord_out = StandardizeCoord(cfg,Coord_in)
%
% standardizes Coord to bins with a specific size
%
% INPUTS
%  coord_in: raw linearized path (obtained from MakeCoord() or its
%  task-specific variants)
%
% OUTPUTS
%  coord_out: standardized linearized path where each unit is a bin of
%  specified size
%
% config options
%  cfg_def.run_dist = 300; % actual length of Coord (in cm)
%  cfg_def.binsize = 3; % desired bin size (in cm)
%  cfg_def.verbose = 1; % 1 display command window text, 0 don't
%
% MvdM 2015-10-01 initial version

cfg_def.run_dist = 300; % actual length of Coord (in cm)
cfg_def.binsize = 3; % desired bin size (in cm)
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose; fprintf('%s: standardizing %d cm coord to %d cm bin size...\n',mfun,cfg.run_dist,cfg.binsize); end

nBins = round(cfg.run_dist/cfg.binsize);

Coord_out(1,:) = interp1(1:size(Coord_in,2),Coord_in(1,:),linspace(1,size(Coord_in,2),nBins),'linear');
Coord_out(2,:) = interp1(1:size(Coord_in,2),Coord_in(2,:),linspace(1,size(Coord_in,2),nBins),'linear');

