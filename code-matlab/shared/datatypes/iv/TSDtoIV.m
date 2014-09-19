function iv_out = TSDtoIV(cfg_in,tsd_in)
% function iv = TSDtoIV(cfg,tsd)
%
% create interval data from tsd by tresholding
%
% INPUTS:
%
% tsd_in: input tsd
%
% CFG OPTIONS:
% cfg.method = 'zscore';
% cfg.threshold = 5;
% cfg.dcn =  '>'; % '<', '>'
% cfg.merge_thr = 0.05; % merge events closer than this
% cfg.target = []; % which data (label) to use
% cfg.minlen = 0.05; % minimum interval length
%
% OUTPUTS:
%
% iv_out: output interval data
%
% MvdM 2014-06-24

cfg.method = 'zscore';
cfg.threshold = 5;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.target = [];
cfg.minlen = 0.05; % minimum interval length

ProcessConfig; % should take whatever is in cfg_in and put it into cfg!
mfun = mfilename;

iv_out = iv; % initialize new iv struct

% check if conditions are in place
nData = size(tsd_in.data,1);
if nData > 1
    if ~isempty(cfg.target)
        temp_data = getd(tsd_in,cfg.target);
    else
       error('Multiple data dimensions exist but no label is specified.'); 
    end
else
    temp_data = tsd_in.data;
end

% zscore
switch cfg.method
    case 'zscore'
        temp_data = zscore(temp_data);
end

% detect crossings
switch cfg.dcn
    case '>'
        detec = temp_data > cfg.threshold;
    case '<'
        detec = temp_data < cfg.threshold;
end

% pad the detection so we can deal with cases where first or last samples are detects
detec = cat(2,0,detec,0);

dfs = diff(detec);
up_idx = find(dfs == 1);
down_idx = find(dfs == -1) - 1;

up_t = tsd_in.tvec(up_idx);
down_t = tsd_in.tvec(down_idx);

% remove doubles
d = up_t(2:end)-down_t(1:end-1);

merge_idx = find(d < cfg.merge_thr);
down_t(merge_idx) = [];
up_t(merge_idx+1) = [];

% remove too short
iv_len = down_t-up_t;
keep_idx = iv_len > cfg.minlen;

iv_out.tstart = up_t(keep_idx);
iv_out.tend = down_t(keep_idx);

% ensure column vectors
if ~iscolumn(iv_out.tstart)
    iv_out.tstart = iv_out.tstart';
end

if ~iscolumn(iv_out.tend)
    iv_out.tend = iv_out.tend';
end

% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_out.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_out.cfg.history.cfg,{cfg});