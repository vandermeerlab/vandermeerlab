function [iv_out,thr_out] = TSDtoIV(cfg_in,tsd_in)
% function [iv_out,thr_out] = TSDtoIV(cfg,tsd)
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
% cfg.verbose = 1; 1 display command window text, 0 don't
%
% OUTPUTS:
%
% iv_out: output interval data
% thr_out: threshold used (in same units as input; helpful if z-scored
%   threshold is to be applied later to other data)
%
% MvdM 2014-06-24

cfg_def.method = 'zscore';
cfg_def.threshold = 5;
cfg_def.dcn = '>'; % return intervals where threshold is exceeded
cfg_def.merge_thr = 0.05; % merge events closer than this
cfg_def.target = [];
cfg_def.minlen = 0.05; % minimum interval length
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % should take whatever is in cfg_in and put it into cfg!

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

% apply transform to data if requested
thr_out = cfg.threshold;
switch cfg.method
    case 'zscore'
        [temp_data,mu,sigma] = zscore(temp_data);
        switch cfg.dcn
            case '>' 
                thr_out = mu + cfg.threshold.*sigma; % store threshold for output argument
            case '<'
                thr_out = mu - cfg.threshold.*sigma;
        end
    case 'percentile'
        temp_data_sorted = sort(temp_data,'ascend');
        thr_out_idx = round(cfg.threshold.*length(temp_data));
        thr_out = temp_data_sorted(thr_out_idx);
        
        temp_data = tiedrank(temp_data)./length(temp_data); % assign percentile to each data point    
        
end

% detect crossings
switch cfg.dcn
    case '>'
        detec = temp_data > cfg.threshold;
    case '<'
        detec = temp_data < cfg.threshold;
    case 'range'
        detec = temp_data > cfg.threshold(1) & temp_data < cfg.threshold(2);
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

if cfg.verbose
    disp([mfun,': ',num2str(length(iv_out.tstart)),' intervals found.'])
end

% housekeeping
iv_out.cfg.history.mfun = cat(1,iv_out.cfg.history.mfun,mfun);
iv_out.cfg.history.cfg = cat(1,iv_out.cfg.history.cfg,{cfg});