function ts_out = ts(timestamps,varargin)
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
% aacarey edit Dec 2017 fixed bug in nargin == 2, changed nargin == 1 label
%         to char instead of numeric

ts_out.type = 'ts';
ts_out.t = {};
ts_out.label = {};

if nargin == 0
    % Do nothing, output empty TS struct
elseif nargin == 1
    if isnumeric(timestamps) && length(timestamps) == 1
        ts_out.t = cell(1,timestamps);
        ts_out.label = cell(1,timestamps);
        
        for iC = 1:timestamps
            ts_out.t{iC} = zeros(0,1);
            ts_out.label{iC} = num2str(iC);
        end
    elseif iscell(timestamps)
        ts_out.t = timestamps;
        ts_out.label = cell(size(timestamps)); % the ts datatype checker only accepts TS where the .t and .labels are the same length
    else
        error('Input is neither a numeric or cell type.')
    end    

elseif nargin == 2
    % initialization of a ts struct with specified data
    labels = varargin{1};    
    if ~iscell(timestamps)
        error('Input .t data must be in cell format.')
    elseif ~iscell(labels)
        error('Input .label data must be in cell format.')
    end
    ts_out.t = timestamps;
    ts_out.label = labels;
    
else % Too many inputs arguments
    error('Too many input arguments.')
end

% housekeeping
ts_out.cfg.history.mfun{1} = mfilename;
ts_out.cfg.history.cfg{1} = [];