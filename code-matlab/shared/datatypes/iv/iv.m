function iv_out = iv(varargin)
% function iv_out = iv(varargin)
%
% constructor for iv (interval) struct
%
% MvdM 2014-06-24

iv_out.tstart = [];
iv_out.tend = [];

if nargin == 1
    
    if numel(varargin{1}) == 2
       iv_out.tstart = varargin{1}(1);
       iv_out.tend = varargin{1}(2);
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