function iv_out = iv(varargin)
% IV Interval datatype constructor
%   Interval, or iv, data is one of the main datatypes. Each iv struct
%   contains a set of intervals and accompanying information (in usr field).
%   Many functions require or assume that the intervals are non-overlapping.
%
%   function iv_out = IV([tstart tend])
%
%   function iv_out = IV(tstart,tend)
%
%   INPUTS:
%
%     [tstart end]: [nx2] double containing interval start times and end times
%
%          OR
%  
%      tstart: [nx1] double or [1xn] double containing interval start times
%        tend: [nx1] double or [1xn] double containing interval end times 
%
%   OUTPUTS
%      iv_out: iv struct with fields:
%           .type   - 'iv'; datatype identification
%           .tstart - [nx1] double containing interval start times
%           .tend   - [nx1] double containing interval end times
%           .usr    - [] initialized field that can be used to store data
%                     or information pertaining to individual intervals
%                     within the iv struct. For example, signal power can be
%                     stored as iv1.usr.power = [nx1] double (NOTE: usr data
%                     must have the same dimensions as the tstart and tend 
%                     fields).
%           .cfg    - record of the data's history including config 
%                     parameters and functions visited
%
%
% MvdM 2014-06-24
% aacarey edit, Nov 2015

iv_out.type = 'iv';
iv_out.tstart = [];
iv_out.tend = [];
iv_out.usr = [];

if nargin == 1
    
    if size(varargin{1},2) == 2
       iv_out.tstart = varargin{1}(:,1);
       iv_out.tend = varargin{1}(:,2);
    else
       error('Single input argument must have length 2 (tstart, tend)');  
    end
    
elseif nargin == 2
    
    if numel(varargin{1}) == numel(varargin{2})
        iv_out.tstart = varargin{1};
        iv_out.tend = varargin{2};
    else
        error('Input arguments must have same length (tstart, tend)');  
    end
    
end

% ensure column vectors
if ~iscolumn(iv_out.tstart)
    iv_out.tstart = iv_out.tstart';
end

if ~iscolumn(iv_out.tend)
    iv_out.tend = iv_out.tend';
end

% housekeeping
iv_out.cfg.history.mfun{1} = mfilename;
iv_out.cfg.history.cfg{1} = [];