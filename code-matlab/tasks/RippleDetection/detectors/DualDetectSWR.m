function SWR_iv = DualDetectSWR(cfg,CSC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
% aacarey 13 Dec 2017

% Parse cfg parameters
cfg_def.verbose = 1;
cfg_def.rippleband = [140 220];
cfg_def.upper = 3;
cfg_def.lower = 0;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

% Detect SWRs
cfg_temp = [];
cfg_temp.rippleband = cfg.rippleband;
cfg_temp.weightby = 'amplitude';
cfg_temp.smooth = 0; % don't smooth, for hitrate
%cfg_temp.kernel = 'wizard'; % which kernel (if empty, goes to default)
cfg_temp.verbose = 0; % talk to me or not

SWR = OldWizard(cfg_temp,CSC);

% Apply smoothing for lower threshold tstart-tend boundaries
cfg_temp = [];
cfg_temp.where = 'all';
cfg_temp.kernel = 'wizard';
cfg_temp.threshold = 0;
cfg_temp.verbose = 0;
SWR_smooth = ConvTSD(cfg_temp,SWR);

% Apply thresholding
% -- upper
cfg_temp =[];
cfg_temp.method = 'zscore';
cfg_temp.threshold = cfg.upper;
cfg_temp.operation =  '>'; % return intervals where threshold is exceeded
cfg_temp.ResizeAt = [];
cfg_temp.target = [];
cfg_temp.verbose = 0;
SWR_detected = TSDtoIV2(cfg_temp,SWR);

% -- lower
cfg_temp.threshold = cfg.lower;
SWR_boundaries = TSDtoIV2(cfg_temp,SWR_smooth);

% Select intervals in SWR_boundaries that are represented in SWR_detected
cfg_temp = [];
cfg_temp.verbose = 0;
SWR_iv = OverlapIV(cfg_temp,SWR_boundaries,SWR_detected);

end

