function TSD_out = ConvTSD(cfg,TSD)
%CONVTSD Convolve time stamped data
%
%   TSD = CONVTSD(cfg,TSD) reshapes TSD.data by convolution with a selected
%   kernel, with the option to target everything or only peaks or valleys. 
%   The intended purpose of this function is to provide a way to emphasize
%   certain features in TSD, especially SWR and MUA detectors.
%   For example, using a sharpness kernel can help accentuate boundaries
%   in SWR detectors, while targeting valleys can help prevent merging of
%   nearby ripples when thresholding to produce IV.
%
%    INPUTS
%      cfg: config struct with fields controlling function behaviour (see
%            CONFIG OPTIONS)
%      TSD: time stamped data to be convolved (.data field only)
%
%    OUTPUTS:
%      TSD: convolved time-stamped data
%
%    CONFIG OPTIONS
%      cfg.where = 'all'; Which regions of TSD.data are being convolved
%           'all'     - Convolve entire.
%           'peaks'   - Convolve peaks only. Kernels like 'wizard' or
%                       'fang' are recommended for ripples.
%           'valleys' - Convolve valleys only. Regions where valleys were
%                       detected are set to zero. Narrow kernels like
%                       'thorn' or 'needle' are recommended for ripple 
%                       separation.
%
%      cfg.kernel = 'wizard'; Kernel that will be convolved with the data
%           'wizard'       - preset gausskernel(60,20)
%           'fang'         - preset gausskernel(200,10)
%           'thorn'        - preset gausskernel(20,2)
%           'needle'       - preset gausskernel(10,1)
%           'triangle'     - preset triangular kernel
%           'sharp1'       - preset sharpness kernel
%           'sharp2'       - preset sharpness kernel
%           [n x 1] double - custom smoothing kernel input by user
%
%      cfg.threshold = 0; Threshold, in SDs above mean, for a peak or
%                    valley to be considered for convolution.
%
%      cfg.verbose = 1; If 1, print informative text to the command window;
%                    if 0, be silent.
%
% see also: triclops
%
% aacarey Feb 2016

% Set config defaults
cfg_def.where = 'all';
cfg_def.kernel = 'wizard';
cfg_def.threshold = 0;
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

if ~CheckTSD(TSD); error('TSD is poorly formed, see the tsd constructor function tsd()'); end

if cfg.verbose; fprintf('%s: Reshaping time-stamped data...\n',mfun); end

% handle custom kernel, if applicable
if isnumeric(cfg.kernel)
    kernel = cfg.kernel;
    cfg.kernel = 'custom';
elseif ~isa(cfg.kernel,'char')
    error('cfg.kernel must either be a string specifying a preset kernel, or a [nx1] double specifying a custom kernel')
end

% assign kernel
switch cfg.kernel
    case 'wizard'
        kernel = gausskernel(60,20);
    case 'fang'
        kernel = gausskernel(200,10);
    case 'thorn'
        kernel = gausskernel(20,2);
    case 'needle'
        kernel = gausskernel(10,1);
    case 'triangle'
        kernel = ([1:21,20:-1:1]./441)';
    case 'sharp1' % sharpness detector
        fatness = 30;
        kernel = [0:-1:-fatness,-fatness:2*fatness,2*fatness:-1:-50,-fatness:0];
        %kernel = diff(diff(gausskernel(20,5).*-1));
        kernel = kernel ./ sum(kernel); % normalize to sum to 1
    case 'sharp2'
        %kernel = (-diff(diff(gausskernel(130,35))))+(1/130);
        %kernel = (-diff(diff(gausskernel(70,35))))+(1/70);
        kernel = sum(-diff(diff(gausskernel(60,20))))+(1/60);
        kernel = kernel ./ sum(kernel); % normalize to sum to 1
    case 'custom'
        % was handled already before the switch case
    otherwise
        error('Unrecognized config option specified in cfg.kernel.')
end

% indicator function will be used to reshape the TSD (peaks and valleys)
indicator = tsd(TSD.tvec,zeros(size(TSD.data)));

switch cfg.where
    case 'all'
        % convolve whole thing
        TSD.data = conv(TSD.data,kernel,'same');
        cfg.threshold = []; % not used in this case
        
    case {'peaks','valleys'}
        
        if strcmp(cfg.where,'peaks')
            direction = 1;
            newmin = min(TSD.data); newmax = max(TSD.data);
        else
            direction = -1;
            newmin = max(TSD.data); newmax = min(TSD.data);
        end
        
        % get zscore
        zscore_data = zscore(TSD.data);
        
        % find local extremes
        % when we're doing valleys, direction is -1 which inverts the TSD
        % so we can use findpeaks to find valleys 
        [value,loc] = findpeaks(zscore_data*direction);
        aboveThreshold = abs(value) > cfg.threshold;
        indicator.data(loc(aboveThreshold)) = direction;
        
        % convolve the indicator
        indicator.data = conv(indicator.data,kernel,'same');
        
        % rescale indicator so it is on the same scale as TSD
        indicator.data = rescale(indicator.data,direction*newmin,direction*newmax);
        
    otherwise
        error('Unrecognized config option specified in cfg.where')
end

% combine TSD with the indicator function (this does something only if we
% took the peaks-or-valleys switchcase)
cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.method = 'sum';
TSD_out = MergeTSD(cfg_temp,TSD,indicator);

% if working with valleys, cut off the extreme dips below zero
if strcmp(cfg.where,'valleys')
   TSD_out.data(TSD_out.data < 0) = 0;
end

% give output the same scale as the input (is this sketchy? is this whole function sketchy?)
TSD_out.data = rescale(TSD_out.data,min(TSD.data),max(TSD.data));

% write config history
TSD_out.cfg.history = TSD.cfg.history;
TSD_out = History(TSD_out,mfun,cfg);

end
