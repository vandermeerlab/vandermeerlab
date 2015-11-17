function tsd_out = tsd(varargin)
% TSD Time-stamped data datatype constructor
%   Time-stamped data, or tsd, is one of the main data types. Each tsd
%   struct contains a time vector and corresponding data (such as voltages 
%   or position). 
%
%   function tsd_out = TSD(tvec,data)
%
%   function tsd_out = TSD(tvec,data,'label')
%
%   INPUTS:
%      tvec: time vector (time stamps in ascending order)
%      data: measurements taken at each time stamp
%     label: unique identifier
%
%   OUTPUTS
%      iv_out: iv struct with fields:
%           .type   - 'tsd'; datatype identification
%           .tvec   - [nx1] double, where n is the number of time samples taken
%           .data   - [1xn] double containing measures taken at each sample
%                     point
%           .label  - Unique identifier for the data contained with the tsd,
%                     such as a record of the data source (in the case of a 
%                     CSC, the label corresponds to the tetrode that recorded 
%                     the signal)
%           .cfg    - record of the data's history including config 
%                     parameters and functions visited
%
% MvdM 2014-06-17
% aacarey edit, Nov 2015

tsd_out.type = 'tsd';
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