function tc_out = tc(varargin)
%% TC TUNING CURVES
% 
%   Constructor for tc (tuning curve) struct.      
%   
% UNDER CONSTRUCTION
%
% youkitan 2014-10-27

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
   
end %iterate place cells
   
    
% housekeeping
tc_out.cfg.history.mfun{1} = mfilename;
tc_out.cfg.history.cfg{1} = [];