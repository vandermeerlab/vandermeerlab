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
% MvdM 2016-01-07 faster version based on nearest_idx3 (taken from
%  restrict2())

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
    
    keep = [];
    
    for iT = 1:length(iv_use.tstart)
        
        ind_start = nearest_idx3(iv_use.tstart(iT),in.tvec,1);
        ind_end = nearest_idx3(iv_use.tend(iT),in.tvec,-1);
        keep = [keep ind_start:ind_end];

    end
    keep = unique(keep);
    
elseif isfield(in,'t') % ts
    
    for iC = length(in.t):-1:1
        keep{iC} = [];
            
        lookups = in.t{iC};
        if isempty(lookups)
            continue
        end
                
        for iT = 1:length(iv_use.tstart) %for each input interval...

            iv_start = iv_use.tstart(iT);
            iv_end = iv_use.tend(iT);
              
            ind_start = nearest_idx3(iv_start,lookups,1);
            ind_end = nearest_idx3(iv_end,lookups,-1);
            
            keepvals_new = ind_start:ind_end;

            if (ind_start == 1 && ind_end == 1) ||...
                    (ind_start == length(lookups) && ind_end == length(lookups))
                keepvals_new = find(in.t{iC} >= iv_use.tstart(iT) & in.t{iC} <= iv_use.tend(iT))';
            end
               
            if ind_start <= ind_end
                keep{iC} = [keep{iC} keepvals_new];
            end
            
        end % iterate events
        
    end % iterate cells
    
end

idx = keep;