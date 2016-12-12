function ts_out = ts(varargin)
% TS Timestamp datatype constructor
% The ts datatype differs from tsd in that it contains only the timestamps of
% events with no additional accompanying data.
%
% function ts_out = TS(varargin)
%
%    OUTPUTS:
%      ts_out: [1x1] ts struct with fields:
%           .type   - 'ts'; datatype identification
%           .t      - {1xn} cell, each containing a [nTimestamps x 1] double
%           .label  - {1xn} cell containing unique identifiers for the 
%                     timestamps contained within the .t field, such as a 
%                     record of the source (in the case of S, the label 
%                     corresponds to the tetrode that recorded the spike 
%                     trains)
%           .cfg    - record of the data's history including config 
%                     parameters and functions visited 
%
% see also tsd
%
% MvdM 2014-06-17
% aacarey edit Nov 2015
% youkitan edit 2016-11-22 added initialization with size

ts_out.type = 'ts';
ts_out.t = {};
ts_out.label = {};

if nargin == 1 && isnumeric(varargin{1}) && size(varargin{1},1) == 1 && size(varargin{1},2) == 1
    ts_out.t = cell(1,varargin{1});
    ts_out.label = cell(1,varargin{1});
    
    for iC = 1:varargin{1}
        ts_out.t{iC} = zeros(0,1);
        ts_out.label{iC} = iC;
    end
end

% housekeeping
ts_out.cfg.history.mfun{1} = mfilename;
ts_out.cfg.history.cfg{1} = [];