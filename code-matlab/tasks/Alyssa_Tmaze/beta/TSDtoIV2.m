function iv_out = TSDtoIV2(cfg_in,tsd_in)
% TSDTOIV2 Create interval data from tsd by thresholding.
%
% iv = TSDTOIV2(cfg,tsd)
%
%   INPUTS:
%        cfg: config struct with fields controlling function behavior
%     tsd_in: tsd struct
%
%   CFG OPTIONS:
%
%     cfg.method = 'zscore';
%         'raw' - threshold raw input data 
%      'zscore' - threshold based on standard deviations above the sample
%                 mean
%
%     cfg.threshold = 0; Single number or [1 x 2] double (see cfg.operation)
%                    specifying threshold value(s).
%
%     cfg.operation =  '>'; 
%             '>' - data > threshold
%             '<' - data < threshold
%            '><' - data > threshold(1) & data < threshold(2)
%            '<>' - data < threshold(1) & data > threshold(2)
% 
%         add '=' to any of the above strings:
%        ex: '>=' - data >= threshold
%           '><=' - data >= threshold(1) & data <= threshold(2)
%
%     cfg.ResizeAtMean = 0; If 1, redraw the intervals so that the
%                   boundaries sit at the mean of tsd_in.data (often results
%                   in merging of nearby events); If 0, don't
%                   (this is done for ripple detection in some papers such
%                   as Hippocampal replay of extended experience)
%                   Note: think about whether this makes sense with respect
%                   to your other config options.
%     
%     cfg.target = []; Which data (label) to use
%     cfg.verbose = 1; 1 - tell me how many intervals you found, 0 - don't
%
%   OUTPUT:
%     iv_out: output interval data
%
%  see also iv TSDtoIV 
%
% MvdM 2014-06-24
% aacarey edit Oct 2015, Feb 2016

% set cfg defaults
cfg_def.method = 'zscore';
cfg_def.threshold = 0;
cfg_def.operation =  '>'; % return intervals where threshold is exceeded
cfg_def.ResizeAtMean = 0;
cfg_def.target = [];
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); 

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
    case 'raw'
        % do nothing
    otherwise
        error('Unrecognized cfg.method.')
end

% detect crossings
switch cfg.operation
    case '>'
        detec = temp_data > cfg.threshold;
    case '>='
        detec = temp_data >= cfg.threshold;
    case '<'
        detec = temp_data < cfg.threshold;
    case '<='
        detec = temp_data <= cfg.threshold;
        
    % the following are kind of unecessary since you can just do any
    % combination of the above in separate calls to TSDtoIV and then put
    % them in a single iv struct.
    case '><'
        detec = temp_data > cfg.threshold(1) & temp_data < cfg.threshold(2);
    case '><='
        detec = temp_data >= cfg.threshold(1) & temp_data <= cfg.threshold(2);
    case '<>'
        detec = temp_data < cfg.threshold(1) & temp_data > cfg.threshold(2);
    case '<>='
        detec = temp_data <= cfg.threshold(1) & temp_data >= cfg.threshold(2);
    otherwise
        error('Unrecognized cfg.operation')
end

% pad the detection so we can deal with cases where first or last samples are detects
detec = [0,detec,0];

dfs = diff(detec);
up_keep = dfs == 1;
down_keep = dfs == -1;

% account for padding (added extra samples)
down_keep = circshift(down_keep,[0,-1]); down_keep(end) = 0;

% make output
up_t = tsd_in.tvec(up_keep);
down_t = tsd_in.tvec(down_keep);

iv_out = iv(up_t,down_t);

if cfg.ResizeAtMean
    % this implementation is probably slow, because the logical operations
    % in OverlapIV are slow.
    
    % threshold at the mean
    cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.method = 'zscore'; cfg_temp.threshold = 0; cfg_temp.operation = '>=';
    iv_mean = TSDtoIV2(cfg_temp,tsd_in); % so meta
    
    % Keep intervals in iv_mean that overlap with intervals in iv_out
    cfg_temp = []; cfg_temp.verbose = 0;
    iv_out = OverlapIV(cfg_temp,iv_mean,iv_out);
   
end

if cfg.verbose
    disp([mfun,': ',num2str(length(iv_out.tstart)),' intervals found.'])
end

% write history
iv_out.cfg.history = tsd_in.cfg.history;
iv_out = History(iv_out,mfun,cfg);
