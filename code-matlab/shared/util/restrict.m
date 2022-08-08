function in = restrict(in,varargin)
% function out = restrict(in,varargin)
%
% restricts times in data object to specific intervals
%
% if a data time is equal to an interval edge, this data point is INCLUDED
%
% usage:
% out = restrict(in,iv)
% out = restrict(in,tstart,tend) % tstart and tend can be vectors
%
% works on ts, tsd, iv data
%
% note that for iv data, intervals must fit within restrict times (i.e.
% will not get cut)
%
% MvdM 2014-07-20 initial version
% youkitan 2016-11-27 edit: type checking update

% convert input arguments to iv if not already done
if nargin == 2
    if ~strcmp(varargin{1}.type,'iv')
       error('Single input argument must be iv (interval) type.'); 
    end
    
    iv_use = varargin{1};
    
elseif nargin == 3
    
    iv_use = iv(varargin{1},varargin{2});
    
else
   error('Unsupported number of input arguments.'); 
end
 
% get indices to keep
if isfield(in,'tvec') && ~isfield(in,'tstart') % tsd
    type = 'tsd';
    keep = false(size(in.tvec));
    
    for iT = 1:length(iv_use.tstart)
        keep = keep | (in.tvec >= iv_use.tstart(iT) & in.tvec <= iv_use.tend(iT));
    end
    
elseif isfield(in,'tstart') % iv
    type = 'iv';
    keep = false(size(in.tstart));
    
    for iT = 1:length(iv_use.tstart)
        keep = keep | (in.tstart >= iv_use.tstart(iT) & in.tend <= iv_use.tend(iT));
    end
    
elseif isfield(in,'t') % ts
    type = 'ts';
    for iC = length(in.t):-1:1
        keep{iC} = false(size(in.t{iC}));
        
        for iT = 1:length(iv_use.tstart)
            keep{iC} = keep{iC} | (in.t{iC} >= iv_use.tstart(iT) & in.t{iC} <= iv_use.tend(iT));
        end
    end
    
end

% do the right thing depending on data type
switch type
    
    case 'ts'
        for iC = 1:length(in.t)
            in.t{iC} = in.t{iC}(keep{iC});
            
            % fixes issues with weird 0x0 cells
            if isempty(in.t{iC})
                in.t{iC} = zeros(0,1);
            end
        end
        
    case 'tsd'
        in.tvec = in.tvec(keep);
        in.data = in.data(:,keep);
    
    case 'iv'
        in.tstart = in.tstart(keep);
        in.tend = in.tend(keep);
        
end

% housekeeping
cfg = []; cfg.iv = iv_use;
in.cfg.history.mfun = cat(1,in.cfg.history.mfun,mfilename);
in.cfg.history.cfg = cat(1,in.cfg.history.cfg,{cfg});