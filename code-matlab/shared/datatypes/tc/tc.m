function tc_out = tc(varargin)
%% TC TUNING CURVE DATATYPE
%   The tc dataype is a container for tuning curves typically used in experimental
%   analysis.
%
%   tc_out = TC(varargin) is a constructor for the tc (tuning curve) struct.
%
%   INPUTS:
%       tc/tc2D: firing rate map for each cell (spike count normalized by sampling of
%       tuning variable)
%       occ_hist: total sampling of tuning variable (occupancy for
%       position variable) spk_hist: binned raw firing rates (optionally smoothed)
%
%   OUTPUTS
%       tc_out: tc struct with fields:
%          .tc - [NxM] array, where N is the number of cells and M is the number of
%                tuning variable bins. 
%        .tc2D - [NxMxP] array, where N is the number of cells and M and P are the number
%                of tuning variable bins(for 2D tuning_var only)
%    .occ_hist - [1xM] array, where M is the number of tuning variable bins
%    .spk_hist - [1xM] array, where M is the number of tuning variable bins
%   
%   See also, TuningCurves, MakeTC, DetectPlaceCells1D
%   Workflow example WORKFLOW_PlotOrderedRaster
%
% youkitan init 2014-10-27
% youkitan edit 2015-09-16
% youkitan edit Feb 2017 overhaul

%% Set fields
tc_out.tc = []; %nCells x nBins
tc_out.occ_hist = []; %
tc_out.spk_hist = [];

switch nargin
    case 2
        tc_out.tc = varargin{1};
        tc_out.occ_hist = varargin{2};
    
    case 3
        tc_out.tc = varargin{1};
        tc_out.occ_hist = varargin{2};
        tc_out.spk_hist = varargin{3};

    case 4
        tc_out.tc = varargin{1};
        tc_out.occ_hist = varargin{2};
        tc_out.spk_hist = varargin{3};
        tc_out.usr.binEdges = varargin{4};
end

% housekeeping
tc_out.cfg.history.mfun{1} = mfilename;
tc_out.cfg.history.cfg{1} = [];