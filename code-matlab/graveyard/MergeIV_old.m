function iv_out = MergeIV(cfg_in,iv1,iv2,varargin)
% iv = MergeIV(cfg,iv1,iv2) returns the union of all input ivs with their overlapping
% portions combined. If intersectFlag is true, returns the true intersect of all input ivs.
%
%   INPUT:
%       cfg_in - input cfg parameters
%       iv1,2 - input ivs (minimum of two ivs)
%       varargin - for more than 2 ivs
%
%   OUTPUT:
%       iv_out - iv containing merged input ivs
%
%   CFG OPTIONS:
%       cfg.dt = 0.05;
%       cfg.mergebin_thr = 5; % in bins, for merge overlaps
%       cfg.mergeFlag = 0; % merge overlaps
%       cfg.intersectFlag = 0; %use intersects of ivs
%
% youkitan 2014-12-31

%% Parse cfg parameters
cfg_def = [];
cfg_def.dt = 0.05;
cfg_def.mergebin_thr = 5;
cfg_def.mergeFlag = 0;
cfg_def.intersectFlag = 0; %use intersects of ivs

cfg = ProcessConfig2(cfg_def,cfg_in);
mfun = mfilename;

%% Make one big iv using UnionIV

% combine iv1 and iv2
cfg_temp = [];
stack_iv = UnionIV(cfg_temp,iv1,iv2);
numIV = 2;

%Check varargin and combine with stack
if nargin > 3
    for iV = 1:length(varargin)
        %check to make sure input is iv
        assert(isfield(varargin{iV},'tstart'),'input must be of iv datatypes') 
    
        %combine iv with stack
        stack_iv = UnionIV([],stack_iv,varargin{iV});
        
        % add to number of IVs
        numIV = numIV + 1;
    end %iterate additional ivs
end

%% Merge overlapping time intervals
lenIV = length(stack_iv.tstart);

% Create internal tsd
tvec_start = stack_iv.tstart(1);
tvec_end = stack_iv.tend(end);
temp_tvec = tvec_start:cfg.dt:tvec_end;
temp_tsd = tsd;
temp_tsd.tvec = temp_tvec;
temp_tsd.data = zeros(1,length(temp_tvec));

% Project iv time intervals onto internal tsd
for iIV = 1:lenIV
    curr_tstart = stack_iv.tstart(iIV);
    curr_tend = stack_iv.tend(iIV);
    
    keep_idx = find(curr_tstart <= temp_tsd.tvec & temp_tsd.tvec <= curr_tend);
    temp_tsd.data(keep_idx) = temp_tsd.data(keep_idx)+1; %timebins that align with ivs
end

if cfg.mergeFlag
    % Merge nearby time intervals (currently super messy)
    temp_ivs = find(temp_tsd.data > 0);
    temp_diff = diff(temp_ivs);
    merge_idx = find(temp_diff > 1 & temp_diff < cfg.mergebin_thr);
    keepiv_idx = [];
    for iM = 1:length(merge_idx)
        numidx = temp_diff(merge_idx(iM));
        temp_start = temp_ivs(merge_idx(iM));
        for iN = 1:numidx
            keepiv_idx = cat(1,keepiv_idx,temp_start+iN);
        end
    end
    temp_tsd.data(keepiv_idx) = 1;
end
    
% Create new iv from tsd
cfg_temp = [];
cfg_temp.method = 'raw';
cfg_temp.dcn =  '>';

if cfg.intersectFlag
    cfg_temp.threshold = numIV-0.01;
else
    cfg_temp.threshold = 0.99;
end

iv_out = TSDtoIV(cfg_temp,temp_tsd);

% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_out.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_out.cfg.history.cfg,{cfg});