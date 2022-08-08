function [linpos_tsd,edges] = LinearizePos(cfg_in,pos_tsd,Coord_in,varargin)
%% LINEARIZEPOS Linearize 2-D position data by projecting onto an alternate coordinate set
%
%   linpos_tsd = LinearizePos(cfg_in,pos_tsd,Coord_in) linearizes (x,y) position
%   data onto Coord_in
%
%	INPUTS:
%       pos_tsd: input position tsd (must have x,y fields)
%       coord_in: raw linearized path (obtained from MakeCoord() or its task-specific
%       variants)
%
%	OUTPUTS:
%       linpos_tsd: tsd with linearized position (z field) and distance (z_dist)
%
%   CFG OPTIONS:
%       cfg_def.debugMode = 0; calculate distance from pos_tsd to coord
%       cfg_def.verbose = 1;
%       cfg_def.nanMethod = 1; % default 1, use nans; if 2, remove nans
%       cfg_def.outputType = 'idx'; If 'idx', returns index of coord points. If 'dist',
%       returns distance in units of run_dist.
%
%   LINEARIZATION EXAMPLE:
%     a          b
%               c
%    1   2   3   4   5
%         d
%
%   In the above schematic, the letters represent the data in pos_tsd and the numbers
%   represent the data in coord. Each letter is projected onto the line defined by the
%   numbers such that the number closest to each letter can be defined as the index of
%   coord points [1, 2, 4, 4]. This is the output data if outputType = 'idx'. If
%   outputType = 'dist', then it is the relative distance from 0 to run_dist. For example,
%   if run_dist = 10cm then the maximum distance (point 5) is 20cm. So, the possible
%   distances are linspace(0,20,5) = [0,5,10,15,20]. The output data is then [0,5,15,15].
%
%   See also, MakeCoord, StandardizeCoord
%   Workflow example Loading_and_Linearizing_Position_Data
%
% Mvdm init
% youkitan edit Feb 2017, overhaul

%% Parse cfg and error check inputs

cfg_def.debugMode = 0;
cfg_def.verbose = 1;
cfg_def.nanMethod = 1;
cfg_def.outputType = 'idx'; %alternatively, 'dist'

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose; 
    disp('LinearizePos: linearizing data...'); 
end

% check for units
if ~strcmp(Coord_in.units,pos_tsd.units)
    fprintf('%s: Inputs do not have same units. Converting variables to match pos_tsd...\n',mfun)
    
    if strcmp(Coord_in.units,'cm') && strcmp(pos_tsd.units,'px')
        Coord_in.coord(1,:) = Coord_in.coord(1,:).*pos_tsd.cfg.ExpKeys.convFact(1); 
        Coord_in.coord(2,:) = Coord_in.coord(2,:).*pos_tsd.cfg.ExpKeys.convFact(2); 
    elseif strcmp(Coord_in.units,'px') && strcmp(pos_tsd.units,'cm')
        Coord_in.coord(1,:) = Coord_in.coord(1,:)./pos_tsd.cfg.ExpKeys.convFact(1); 
        Coord_in.coord(2,:) = Coord_in.coord(2,:)./pos_tsd.cfg.ExpKeys.convFact(2); 
    else
        error('Currently unsupported unit type.')
    end
    
end

% check for standardization
if ~Coord_in.standardized    
    fprintf('%s: WARNING, input coords are not standardized.\n',mfun);
end

% check for 'dist' output requirements
if strcmp(cfg.outputType,'dist')
    if Coord_in.standardized
        run_dist = Coord_in.run_dist;
    elseif nargin > 3 && strcmp(varargin{1},'run_dist')
        run_dist = varargin{2};
    else
        error('Selected output type requires run_dist input.')
    end
    
    % type check
    if ~isnumeric(run_dist)
        error('Input run_dist should be numeric.')
    end
   
end

%% main body

% init output variable
linpos_tsd = tsd; 
linpos_tsd.units = cfg.outputType;
linpos_tsd.cfg = pos_tsd.cfg;

% set up vars for linearize
x = getd(pos_tsd,'x'); 
y = getd(pos_tsd,'y');

if any(isnan(x) | isnan(y))
    fprintf('%s: NaNs detected in input position. NaN method %d has been selected...\n',mfun,cfg.nanMethod);
end

linpos_temp = NaN(size(x));
keep_idx = ~isnan(x) & ~isnan(y);

% ensures loops don't have duplications
if Coord_in.coord(1,1) == Coord_in.coord(1,end) && Coord_in.coord(2,1) == Coord_in.coord(2,end)
    Coord_in.coord = Coord_in.coord(:,1:end-1);
end

% get Coord points nearest to each x,y pair
coord_vals = 1:length(Coord_in.coord); % set of possible values for linpos_tsd.data
linpos_temp(keep_idx) = griddata(Coord_in.coord(1,:),Coord_in.coord(2,:),coord_vals,x(keep_idx),y(keep_idx),'nearest');

% also get distance
dist_temp = NaN(size(linpos_temp));
dist_temp(keep_idx) = sqrt((Coord_in.coord(1,ceil(linpos_temp(keep_idx))) - x(keep_idx)).^2 + (Coord_in.coord(2,ceil(linpos_temp(keep_idx))) - y(keep_idx)).^2);

% take care of Nans
if cfg.nanMethod == 2 % remove any data points that contain nans
    nan_idx = isnan(x) | isnan(y);
    linpos_temp(nan_idx) = [];
    dist_temp(nan_idx) = [];
    pos_tsd.tvec(nan_idx) = [];
end

switch cfg.outputType
    case 'dist'
        dist_vals = linspace(0,run_dist,length(Coord_in.coord));
        linpos_temp = dist_vals(linpos_temp);
        
        jump = abs(dist_vals(1)-dist_vals(2));
        edges = [0,(dist_vals(2:end) - (jump./2)),(run_dist)]; %first and last bin are a bit weird?
        
    case 'idx'
        edges = min(linpos_temp)-0.5:1:max(linpos_temp)+0.5;
end

% assemble output tsd
linpos_tsd.tvec = pos_tsd.tvec;
linpos_tsd.data(1,:) = linpos_temp;
linpos_tsd.label{1} = 'z';

if cfg.debugMode
    linpos_tsd.data(2,:) = dist_temp;
    linpos_tsd.label{2} = 'z_dist';
end

%% Housekeeping
linpos_tsd = History(linpos_tsd,mfun,cfg);

