function ts_out = ts(varargin)
% TS Timestamp datatype constructor
%   The ts datatype differs from tsd in that it contains only the timestamps of events
%   with no additional accompanying data.
%   
%   function ts_out = TS(t,label)
%
%   function ts_out = TS(t)
%
%   function ts_out = TS(N) returns an ts struct with N specified empty cells
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
% see also tsd, CheckTS
%
% MvdM 2014-06-17
% aacarey edit Nov 2015
% youkitan edit 2016-11-22 added initialization with size
% youkitan edit Jan 2017 added intitialization with data

ts_out.type = 'ts';
ts_out.t = {};
ts_out.label = {};

% initialization of a ts struct with specified data
if nargin == 2
    timestamps = varargin{1};
    labels = varargin{2};
    
    if ~iscell(timestamps)
        error('Input .t data are not in cell format.')
    elseif ~iscell(labels)
        error('Input .label data are not in cell format.')
    end
    
    ts_out.t = timestamps;
    ts_out.label = labels;

% initialization of an empty ts struct with prespecified number of data bins
elseif nargin == 1 
    if isnumeric(varargin{1}) && size(varargin{1},1) == 1 && size(varargin{1},2) == 1
        ts_out.t = cell(1,varargin{1});
        ts_out.label = cell(1,varargin{1});

        for iC = 1:varargin{1}
            ts_out.t{iC} = zeros(0,1);
            ts_out.label{iC} = iC;
        end
    elseif iscell(varargin{1})
        ts_out.t = varargin{1};
    else
        error('Input is neither a numeric or cell type.')
    end

% initialization of an empty ts struct
elseif nargin == 0
    %doesn't have to do aything
else
    error('Invalid number of arguments.')
end

% housekeeping
ts_out.cfg.history.mfun{1} = mfilename;
ts_out.cfg.history.cfg{1} = [];