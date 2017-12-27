function [threshold,tstartThreshold,tendThreshold] = GetThreshold(cfg_in,iv_in,tsd_in)
%GETTHRESHOLD Returns the threshold that will generate intervals like iv_in
%from the data in tsd_in.
% [threshold,tstartThreshold,tendThreshold] = GetThreshold(cfg,iv,tsd)
%
%   INPUTS:
%        cfg: config struct with fields controlling function behavior
%        tsd: tsd struct
%
%   CFG OPTIONS:
%
%     cfg.method = 'zscore';
%         'raw' - uses raw input data to estimate the threshold
%      'zscore' - uses zscored input data to estimate the threshold
%
%     cfg.return = 'median'; How would you like the threshold to be
%                   calculated
%         'mean' - returns the mean threshold
%       'median' - returns the median threshold
%         'mode' - returns the most-occurring threshold
%
%     cfg.verbose = 1; 
%        If 1, prints text to the command window, if 0 doesn't.
%
% aacarey Nov 2017, initial version

%Set cfg defaults
cfg_def.method = 'zscore'; % 'zscore' or 'raw'
cfg_def.return = 'median'; % 'mean', 'median','mode'
cfg_def.verbose = 1;

% Parse cfg parameters
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); 

tvec = tsd_in.tvec;
data = tsd_in.data;

switch cfg.method
    case 'zscore'
        % get zscore
        data = zscore(data);
    case 'raw'
        %do nothing
    otherwise
        error('Unrecognized config option specified in cfg.method')
end

% find where the IV boundaries were drawn with respect to tvec
tstart = data(nearest_idx3(iv_in.tstart,tvec));
tend = data(nearest_idx3(iv_in.tend,tvec));

% pull the threshold out of data
switch cfg.return
    case 'mean'
        tstartThreshold = mean(tstart);
        tendThreshold = mean(tend);
        threshold = mean([tstart tend]);
    case 'median'
        tstartThreshold = median(tstart);
        tendThreshold = median(tend);
        threshold = median([tstart tend]);
    case 'mode'
        tstartThreshold = mode(tstart);
        tendThreshold = mode(tend);
        threshold = mode([tstart tend]);
    otherwise
        error('Unrecognized config option specified in cfg.return')
end

% talk to me or not
if cfg.verbose
   disp([mfun,': the ',cfg.return,' ',cfg.method,' threshold is ',num2str(threshold)]) 
end

end

