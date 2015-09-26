function tc_out = tc(varargin)
%% TC TUNING CURVE DATATYPE
%   tc_out = TC(varargin) is a constructor for the tc (tuning curve) struct.
%
% youkitan 2014-10-27
% youkitan edit 2015-09-16

%% Set fields
tc_out.tc = []; %nCells x nBins
tc_out.occ = []; %
tc_out.label = {};
tc_out.units = {};

switch nargin
    case 2
        tc_out.tc = varargin{1};
        tc_out.tc = varargin{2};
    
    case 3
        tc_out.tc = varargin{1};
        tc_out.tc = varargin{2};
        tc_out.tc = varargin{3};

    case 4
        tc_out.tc = varargin{1};
        tc_out.tc = varargin{2};
        tc_out.tc = varargin{3};
        tc_out.tc = varargin{4};
end

% Enforce equal sized arrays
for i = 1:length(tc_out.tc)
    % Check to make sure vectors are of same size
    assert(length(tc_out.firing_rates{i}) == length(tc_out.stimulus{1,i}),...
        'length of vectors must be equal!')    
   
    % Check to see how many dimensions you have for the stimulus
   dims = size(tc_out.stimulus{i});
   fprintf('The stimulus for cell %d has %d number of dimensions',i,dims(1));
   
end %iterate cells



%% Set defaults
tc_out.firing_rates = {}; %nCells x nBins
tc_out.stimulus = {}; %nCells x nStimulus
tc_out.parameters = [];
tc_out.label = {};

% Need to enforce same stimulus parameters for all cells

switch nargin
    case 2
        tc_out.firing_rates = varargin{1};
        tc_out.stimulus = varargin{2};

    case 3
        tc_out.firing_rates = varargin{1};
        tc_out.stimulus = varargin{2};
        tc_out.label = varargin{3};

    case 4
        tc_out.firing_rates = varargin{1};
        tc_out.stimulus = varargin{2};
        tc_out.parameters = varargin{3};
        tc_out.label = varargin{4};
end

for i = 1:length(tc_out.stimulus)
    % Check to make sure vectors are of same size
    assert(length(tc_out.firing_rates{i}) == length(tc_out.stimulus{1,i}),...
        'length of vectors must be equal!')    
   
    % Check to see how many dimensions you have for the stimulus
   dims = size(tc_out.stimulus{i});
   fprintf('The stimulus for cell %d has %d number of dimensions',i,dims(1));
   
end %iterate cells
   
    
% housekeeping
tc_out.cfg.history.mfun{1} = mfilename;
tc_out.cfg.history.cfg{1} = [];