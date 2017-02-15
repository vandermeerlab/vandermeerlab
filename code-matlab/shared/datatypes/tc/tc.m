function tc_out = tc(varargin)
%% TC TUNING CURVE DATATYPE
%   The tc dataype is a container for tuning curves typically used in experimental
%   analysis.
%
%   tc_out = TC(varargin) is a constructor for the tc (tuning curve) struct.
%
%   INPUTS:
%
%
%   OUTPUTS
%       iv_out: iv struct with fields:
%          .tc - [NxM] array, where N is the number of cells and M is the number of
%                tuning variable bins. Each entry is a normalized firing rate.
%        .tc2D - [NxMxP] array, 2-dimensional firing rate map (for 2D tuning_var only)
%    .occ_hist - [1xM] array, occupancy (sample counts of tuning variable)
%    .spk_hist - [1xM] array, raw spike counts
%
% youkitan init 2014-10-27
% youkitan edit 2015-09-16
% youkitan edit Feb 2017 overhaul

%% Set fields
tc_out.tc = []; %nCells x nBins
tc_out.occ_hist = []; %
tc_out.spk_hist = [];
tc_out.binEdges = [];
tc_out.binCenters = [];
tc_out.nBins = [];

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
        tc_out.tc = varargin{4};
end

% housekeeping
tc_out.cfg.history.mfun{1} = mfilename;
tc_out.cfg.history.cfg{1} = [];