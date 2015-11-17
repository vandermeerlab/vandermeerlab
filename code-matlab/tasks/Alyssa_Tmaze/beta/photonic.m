function tsd_out = photonic(cfg_in,tsd_in)
%PHOTONIC Detect transients in the local field potential
%
%    Transient detection works best with little background noise. Since the 
%  HC has a lot of high-amplitude, low frequency background "noise", you 
%  should highpass or bandpass filter the signal before sending it in (see 
%  FilterLFP and fftFilterLFP).
%
% tsd_out = PHOTONIC(cfg,CSCf)
%
%    INPUTS
%          cfg: config struct with fields controlling function behavior
%       tsd_in: tsd struct
%
%    OUTPUTS
%      tsd_out: tsd struct with .data field containing detection results
%
%    CONFIG OPTIONS 
%       cfg.weightby = 'power'; How to obtain the signal envelope
%          'amplitude'    - amplitude envelope (absolute value of the signal)
%          'power'        - power envelope (square the signal)
% 
%       cfg.kernel = 'gauss'; Which smoothing kernel to use
%          'gauss'        - preset gaussian kernel
%          'triangle'     - preset triangular kernel
%          'sharp1'       - preset sharpness detector
%          'sharp2'       - preset sharpness detector
%          [n x 1] double - custom smoothing kernel
%
%       cfg.verbose = 1; Tell me that you are doing something.
%
% (proposed mundane name for photonic: DetectTransients) 
%  aacarey Oct 2015
% 
% see also OldWizard

%%
mfun = mfilename;

% set cfg defaults
cfg_def.weightby = 'power';
cfg_def.kernel = 'gauss';
cfg_def.verbose = 1;

% parse cfg parameters
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% if ~CheckTSD(lfp_tsd)
%    error('Input must be tsd struct with .data and .tvec fields. See tsd() or LoadCSC().'); 
% end

% talk to me or not
if cfg.verbose
   disp([mfun,': looking for transients in the signal...']) 
end

% obtain the signal's envelope
switch cfg.weightby
    case 'amplitude' % amplitude envelope
        envelope = abs(tsd_in.data);
    case 'power' % power envelope
        envelope = tsd_in.data.^2;
    otherwise
        error('Unrecognized cfg.weightby')
end

% handle custom kernel, if applicable
if isnumeric(cfg.kernel)
    kernel = cfg.kernel;
    cfg.kernel = 'custom';
elseif ~isa(cfg.kernel,'char')
    error('cfg.kernel must either be a string specifying a preset kernel, or a [nx1] double specifying a custom kernel')
end

% assign kernel
switch cfg.kernel
    case 'gauss'
        kernel = gausskernel(60,20);
    case 'triangle'
        kernel = ([1:21,20:-1:1]./441)';
    case 'sharp1' % sharpness detector
        fatness = 30;
        kernel = [0:-1:-fatness,-fatness:2*fatness,2*fatness:-1:-50,-fatness:0];
%         kernel = diff(diff(gausskernel(20,5).*-1));
        kernel = kernel ./ sum(kernel); % normalize to sum to 1
    case 'sharp2'
        %kernel = (-diff(diff(gausskernel(130,35))))+(1/130);
        %kernel = (-diff(diff(gausskernel(70,35))))+(1/70);
        kernel = (-diff(diff(gausskernel(60,20))))+(1/60);
    case 'custom'
        % do nothing because already done ^^
    otherwise
        error('Unrecognized cfg.kernel.')
end

% apply smoothing
smooth = conv(envelope,kernel,'same');

% sharpen
% kernel = diff(diff(gausskernel(100,20).*-1));
% kernel = kernel ./ sum(kernel); % normalize to sum to 1
% smooth = conv(smooth,kernel,'same');

% construct tsd output
tsd_out = tsd(tsd_in.tvec,smooth);

% housekeeping
tsd_out.cfg.history.mfun = cat(1,tsd_in.cfg.history.mfun,{mfun});
tsd_out.cfg.history.cfg = cat(1,tsd_in.cfg.history.cfg,{cfg});

end

