function linpos_tsd = LinearizePos(cfg_in,pos_tsd)
% function linpos_tsd = LinearizePos(cfg,pos_tsd)
%
% linearizes (x,y) position data
%
% INPUTS:
%
% pos_tsd: input position tsd (must have x,y fields)
%
% OUTPUTS:
%
% linpos_tsd: tsd with lineariezed position (z field) and distance (z_dist)
%
% cfg.Coord: Coord file geneated by makeCoord, if not defined tries to load

cfg.Coord = [];
ProcessConfig; % remember this will overwrite Coord if defined in cfg

mfun = mfilename;

disp('LinearizePos.m: linearizing data...');

% try to load Coord file if not specified
if isempty(cfg.Coord)
   
    fn = FindFiles('*Coord.mat','CheckSubdirs',0);
    
    if length(fn) > 1
       error('More than one Coord file found. Specifiy file to use in cfg.'); 
    else
        fn = fn{1};
    end
    
    load(fn);
    cfg.Coord = Coord;
    
end

% init output variable
linpos_tsd = tsd; linpos_tsd.cfg = pos_tsd.cfg;

% set up vars for linearize
x = getd(pos_tsd,'x'); y = getd(pos_tsd,'y');

NN = repmat(NaN,size(x));
f_use = find(~isnan(x) & ~isnan(y));

if cfg.Coord(1,1) == cfg.Coord(1,end) & cfg.Coord(2,1) == cfg.Coord(2,end)
    cfg.Coord = cfg.Coord(:,1:end-1);
end

% get Coord points nearest to each x,y pair
NN(f_use) = griddata(cfg.Coord(1,:),cfg.Coord(2,:),1:length(cfg.Coord),x(f_use),y(f_use),'nearest');

% also get distance
d = sqrt((cfg.Coord(1,ceil(NN(f_use))) - x(f_use)).^2 + (cfg.Coord(2,ceil(NN(f_use))) - y(f_use)).^2);

% assemble output tsd
linpos_tsd.tvec = pos_tsd.tvec;

linpos_tsd.data(1,:) = NN;
linpos_tsd.label{1} = 'z';

linpos_tsd.data(2,:) = d;
linpos_tsd.label{2} = 'z_dist';


% housekeeping
linpos_tsd.cfg.history.mfun = cat(1,linpos_tsd.cfg.history.mfun,mfun);
linpos_tsd.cfg.history.cfg = cat(1,linpos_tsd.cfg.history.cfg,{cfg});

