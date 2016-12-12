function in = restrict2(in,varargin)
%RESTRICT2 function out = restrict2(in,varargin)
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
% youkitan 2015-07-11 version 2
% aacarey edit, Nov 25 2015
% youkitan edit 2016-03-28

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
if isfield(in,'tstart') % iv
    type = 'iv';
    keep = [];
    
    for iT = 1:length(iv_use.tstart)
%         keep = keep | (in.tstart >= iv_use.tstart(iT) & in.tend <= iv_use.tend(iT));
        ind_start = nearest_idx3(iv_use.tstart(iT),in.tstart,1);
        ind_end = nearest_idx3(iv_use.tend(iT),in.tstart,-1);
        keep = [keep ind_start:ind_end];
        
        assert(ind_start <= ind_end,'Indices incorrect!')

    end
    
elseif isfield(in,'tvec') % tsd
    type = 'tsd';
    
    idx_start = nearest_idx3(iv_use.tstart,in.tvec,1);
    idx_end = nearest_idx3(iv_use.tend,in.tvec,-1);
    keep =[];
    
    for iT = 1:length(idx_start)
        keep = [keep idx_start(iT):idx_end(iT)];
    end
    
    % there are multiple copies of some indices, so remove them
    keep = unique(keep);
    
elseif isfield(in,'t') % ts
    type = 'ts';
    for iC = length(in.t):-1:1
        keep{iC} = [];
        
        lookups = in.t{iC};
        if isempty(lookups)
            continue
        end
            
        for iT = 1:length(iv_use.tstart) %for each input interval...

            iv_start = iv_use.tstart(iT);
            iv_end = iv_use.tend(iT);
              
            % makes valid lookups "within" start/end times
            ind_start = nearest_idx3(iv_start,lookups,1);
            ind_end = nearest_idx3(iv_end,lookups,-1);
            
            % accounts for lookup values that equal a start/end time
            if ind_start >= ind_end && ind_start ~= 1
                if lookups(ind_start-1) == iv_start
                    ind_start = ind_start-1;
                end
            end

            keepvals_new = ind_start:ind_end;
            
            if (ind_start == 1 && ind_end == 1) ||...
                    (ind_start == length(lookups) && ind_end == length(lookups))
                keepvals_new = find(in.t{iC} >= iv_use.tstart(iT) & in.t{iC} <= iv_use.tend(iT))';
            end
               
            if ind_start <= ind_end
                keep{iC} = [keep{iC} keepvals_new];
            end
            
        end %iterate events
        
    end %iterate cells
    
end

% do the right thing depending on data type
switch type
    
    case 'ts'
        for iC = 1:length(in.t)
            in.t{iC} = in.t{iC}(keep{iC});
            
            % In original restrict(), empty cells were 0x1 doubles instead of 0x0 doubles.
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