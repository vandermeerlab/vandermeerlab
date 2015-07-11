function in = restrict2(in,varargin)
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
if isfield(in,'tstart') % iv
    
    keep = [];
    
    for iT = 1:length(iv_use.tstart)
%         keep = keep | (in.tstart >= iv_use.tstart(iT) & in.tend <= iv_use.tend(iT));
        ind_start = nearest_idx3(iv_use.tstart(iT),in.tstart,1);
        ind_end = nearest_idx3(iv_use.tend(iT),in.tstart,-1);
        keep = [keep ind_start:ind_end];
        
        assert(ind_start <= ind_end,'Indicies incorrect!')

    end
    
elseif isfield(in,'tvec') % tsd
%     keep = false(size(in.tvec));
    
    keep = [];
    for iT = 1:length(iv_use.tstart)
        %currently takes too long this way O(n^2)
%         keep = keep | (in.tvec >= iv_use.tstart(iT) & in.tvec <= iv_use.tend(iT));
        iv_start = iv_use.tstart(iT);
        iv_end = iv_use.tend(iT);
        lookups = in.tvec;
        
        ind_start = nearest_idx3(iv_start,lookups,-1);
        ind_end = nearest_idx3(iv_end,lookups,-1);

        keepvals_new = [ind_start:ind_end];
%         keepvals_old = find(in.tvec >= iv_use.tstart(iT) & in.tvec <= iv_use.tend(iT))';
%         [not_in_new,not_in_new_ind] = setdiff(keepvals_old,keepvals_new);
%         [not_in_old,not_in_old_ind] = setdiff(keepvals_new,keepvals_old);
        keep = [keep keepvals_new];
        
        assert(ind_start <= ind_end,'Indicies incorrect!')
%         assert(isequal(keepvals_old,keepvals_new),'Different values!\n')
    end
    
elseif isfield(in,'t') % ts
    
    for iC = length(in.t):-1:1
        keep{iC} = [];
        
        for iT = 1:length(iv_use.tstart)
%             keep{iC} = keep{iC} | (in.t{iC} >= iv_use.tstart(iT) & in.t{iC} <= iv_use.tend(iT));

            iv_start = iv_use.tstart(iT);
            iv_end = iv_use.tend(iT);
            lookups = in.t{iC};
            
            ind_start = nearest_idx3(iv_start,lookups,1);
            ind_end = nearest_idx3(iv_end,lookups,-1);
            keepvals_new = [ind_start:ind_end];
            keepvals_old = find(in.t{iC} >= iv_use.tstart(iT) & in.t{iC} <= iv_use.tend(iT))';

            if (ind_start == 1 && ind_end == 1) ||...
                    (ind_start == length(lookups) && ind_end == length(lookups))
                keepvals_new = keepvals_old;
            end
               
            assert(isequal(keepvals_old,keepvals_new),...
                'Different values for cell %d during event %d!\n',iC,iT)
            
            if ind_start <= ind_end
                keep{iC} = [keep{iC} keepvals_new];
            end
            
        end %iterate events
        
    end %iterate cells
    
end

% do the right thing depending on data type
switch in.cfg.history.mfun{1}
    
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