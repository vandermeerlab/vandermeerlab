function [in,keep] = antirestrict(in, varargin)
% function [out,keep] = antirestrict(in,varargin)
%
% cuts out specific intervals from data object
%
% usage:
% out = antirestrict(in,iv)
% out = antirestrict(in,tstart,tend) % tstart and tend can be vectors
%
% see also restrict()
%
% MvdM 2024

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
    
    keep = ~keep;

elseif isfield(in,'tstart') % iv
    type = 'iv';
    keep = false(size(in.tstart));
    
    for iT = 1:length(iv_use.tstart)
        keep = keep | (in.tstart >= iv_use.tstart(iT) & in.tend <= iv_use.tend(iT));
    end

    keep = ~keep;
    
elseif isfield(in,'t') % ts
    type = 'ts';
    for iC = length(in.t):-1:1
        keep{iC} = false(size(in.t{iC}));
        
        for iT = 1:length(iv_use.tstart)
            keep{iC} = keep{iC} | (in.t{iC} >= iv_use.tstart(iT) & in.t{iC} <= iv_use.tend(iT));
        end
        keep{iC} = ~keep{iC};
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
        
        if isfield(in,'usr')
            if ~isempty(in.usr)
               fn = fieldnames(in.usr);
               for iFN = 1:length(fn)
                   temp = in.usr.(fn{iFN});
                   temp = temp(keep);
                   in.usr.(fn{iFN}) = temp;
               end
            end
        end
        
end

% housekeeping
cfg = []; cfg.iv = iv_use;
in.cfg.history.mfun = cat(1,in.cfg.history.mfun,mfilename);
in.cfg.history.cfg = cat(1,in.cfg.history.cfg,{cfg});