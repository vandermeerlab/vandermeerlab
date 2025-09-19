function [in,keep] = restrict(in,varargin)
% function [out,keep] = restrict(in,varargin)
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
% ED 2025-09-19 edit: added possibility to include (cut) first & last IV interval
% when the restrict would have cut it, instead of discarding it. Only
% applies to IV for now.

% see also antirestrict()
keep_ends = 0; % default behaviour: intervals that start before or 
% end after the restricted limits will not be included

% convert input arguments to iv if not already done
if nargin == 2
    if ~strcmp(varargin{1}.type,'iv')
       error('Single input argument must be iv (interval) type.'); 
    end
    
    iv_use = varargin{1};
    
elseif nargin == 3
    
    iv_use = iv(varargin{1},varargin{2});
    
elseif nargin ==4

    iv_use = iv(varargin{1},varargin{2});
    % iv, st, et and "keep_first".
    % if keep_first is 1, will keep the first interval even if it means it
    % would be cut. useful e.g. when dealing with long intervals, like
    % theta mode.
    keep_ends = varargin{3};

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

    % new ED
    if keep_ends
        % Check if there was an interval with end included but start not,
        % cut it to the earliest time
        % find first iv that ends in the interval
        

        % TWO types: 
        % try to find intervals that start before the start and after the start
        % i.e, that would have been removed because they start
        % too early
        early_ivs = find(in.tstart < iv_use.tstart(iT) & in.tend > iv_use.tstart(iT));
        % and those that start before the end, but end after the end
        late_ivs = find(in.tstart < iv_use.tend(iT) & in.tend > iv_use.tend(iT));

        % out of bounds ivs
        oob_ivs = [early_ivs, late_ivs];
        % For each, truncate them
        ivs_to_add = [];
        for iv_i = 1:length(oob_ivs)
            iv_ind = oob_ivs(iv_i);
            % truncate the start to be no earlier than the interval start
            % and the end to be no later than the interval end
            this_st = max(in.tstart(iv_ind), iv_use.tstart(iT));
            this_et = min(in.tend(iv_ind), iv_use.tend(iT));
            ivs_to_add = [ivs_to_add; [this_st, this_et]];
        end        
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


        % new ED
        if keep_ends
            
            % add the truncated intervals to the index to keep
            keep(oob_ivs) = 1;

            % Update the intervals to the truncated version
            for int_i = 1:length(oob_ivs)
                oob_ind = oob_ivs(int_i);
                in.tstart(oob_ind) = ivs_to_add(int_i,1);
                in.tend(oob_ind) = ivs_to_add(int_i, 2);
            end

        end

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