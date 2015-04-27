function idx = TSD_getidx(in,varargin)
% function idx = TSD_getidx(in,varargin)
%
% returns idxs into ts or tsd, delimited by input interval(s)
%
% if a data time is equal to an interval edge, this data point is INCLUDED
%
% usage:
% idx = TSD_getidx(in,iv)
% idx = TSD_getidx(in,tstart,tend) % tstart and tend can be vectors
%
% MvdM 2015-04-26 initial version: SLOW for large numbers of inputs, fix

% convert input arguments to iv if not already done
if nargin == 2
    if ~strcmp(varargin{1}.cfg.history.mfun{1},'iv')
       error('Single input argument must be iv (interval) type.'); 
    end
    
    iv_use = varargin{1};
    
elseif nargin == 3
    
    iv_use = iv(varargin{1},varargin{2});
    
else
   error('Unsupported number of input arguments.'); 
end
 
% get indices to keep
if isfield(in,'tvec') % tsd
    keep = false(size(in.tvec));
    
    for iT = 1:length(iv_use.tstart)
        keep = keep | (in.tvec >= iv_use.tstart(iT) & in.tvec <= iv_use.tend(iT));
    end
    
elseif isfield(in,'t') % ts
    
    for iC = length(in.t):-1:1
        keep{iC} = false(size(in.t{iC}));
        
        for iT = 1:length(iv_use.tstart)
            keep{iC} = keep{iC} | (in.t{iC} >= iv_use.tstart(iT) & in.t{iC} <= iv_use.tend(iT));
        end
    end
    
end

idx = find(keep);