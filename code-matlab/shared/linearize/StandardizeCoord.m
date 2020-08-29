function Coord_out = StandardizeCoord(cfg_in,Coord_in,run_dist,varargin)
%% STANDARDIZECOORD Standardize linearization coordinates
%
%   [Coord_out,cfg_coord] = StandardizeCoord(cfg_in,Coord_in,run_dist) standardizes Coord
%   to 100 points of equal distance by interpolating the raw input Coord.
%
%   [Coord_out,cfg_coord] = StandardizeCoord(cfg_in,Coord_in,run_dist,'nPoints',N)
%   standardizes Coord to N points of equal distance.
%
%   [Coord_out,cfg_coord] = StandardizeCoord(cfg_in,Coord_in,run_dist,'pointDist'',M)
%   standardizes Coord to points equally separated by M units.
%
%   INPUTS
%       coord_in: raw linearized path (obtained from MakeCoord() or its
%       task-specific variants)
%       run_dist: true path length the coord points represent (in experimental units
%
%   OUTPUTS
%       Coord_out - output coord struct
%              .coord: standardized coord data points
%              .units: units of coords (from MakeCoord)
%          . run_dist: input run_dist
%            .nPoints: number of coord points (same as length(.coord))
%          .pointDist: distance between coord points in units of run_dist
%       .standardized: logical flag, is coord standardized or not
%
%   CFG OPTIONS
%       cfg.verbose = 1; If 1, print informative text to the command window; if 0, be
%       silent.
%
% MvdM 2015-10-01 initial version
% youkitan edit Feb 2017, overhaul

%% parse cfg
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose; 
    fprintf('%s: standardizing coord...\n',mfun); 
end

%% main body
Coord_out = Coord_in;

if nargin == 3
    nPoints = 100;
    pointDist = run_dist./nPoints;
    
elseif nargin == 5
    if strcmp(varargin{1},'nPoints')
        nPoints = varargin{2};
        pointDist = run_dist./nPoints;
        
    elseif strcmp(varargin{1},'pointDist')
        pointDist = varargin{2};
        
        % if you specify a distance between points, it needs to factor into the run distance such
        % that nPoints is a whole number
        if mod(run_dist,pointDist) > 0
            fprintf('%s: input pointDist factors incorrectly. Approximating value...\n',mfun)
            nPoints = round(run_dist./pointDist);
            pointDist = run_dist./nPoints;
        else
            nPoints = run_dist./pointDist;
        end
        
    else
        error('Input must be a valid string-value pair. See help.')
    end
else
    error('Input must contain a ParameterName/Value pair. See help.')
end
    
coord_std = [];
coord_std(1,:) = interp1(1:size(Coord_in.coord,2),Coord_in.coord(1,:),linspace(1,size(Coord_in.coord,2),nPoints),'linear');
coord_std(2,:) = interp1(1:size(Coord_in.coord,2),Coord_in.coord(2,:),linspace(1,size(Coord_in.coord,2),nPoints),'linear');

%% hosuekeeping
Coord_out.coord = coord_std;
Coord_out.run_dist = run_dist;
Coord_out.nPoints = nPoints;
Coord_out.pointDist = pointDist;
Coord_out.standardized = 1;

end
