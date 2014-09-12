function tsd_out = tsd(varargin)
% function tsd_out = tsd(varargin)
%
% constructor for tsd (timestamped data) struct
%
% MvdM 2014-06-17

tsd_out.tvec = [];
tsd_out.data = [];
tsd_out.label = {};

if nargin == 2
    tsd_out.tvec = varargin{1};
    tsd_out.data = varargin{2};
end

if nargin == 3
    tsd_out.tvec = varargin{1};
    tsd_out.data = varargin{2};
    tsd_out.label = varargin{3};
end

% housekeeping
tsd_out.cfg.history.mfun{1} = mfilename;
tsd_out.cfg.history.cfg{1} = [];