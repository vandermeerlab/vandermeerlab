function SWR = OldWizard(cfg_in,CSC)
% OLDWIZARD Detect sharp wave-ripples using the Hilbert transform.
%
% SWR = OLDWIZARD(cfg,CSC) filters the signal in the passband specified by
%    cfg.rippleband, then the Hilbert transform is taken and smoothed to
%    create a detection score vector.
%
%   INPUTS:
%         cfg: config struct with fields controlling function behavior
%         CSC:  tsd from a continuously sampled channel (output from LoadCSC) 
%
%   OUTPUTS
%         SWR:  [1 x 1] struct, of tsd datatype, with fields:
%             .cfg - record of cfg history
%             .data - detection score vector [1 x nSamples] double
%             .tvec - time vector from CSC [nTimestamps x 1] double
%
%   CFG OPTIONS
%       cfg.rippleband = [140 250]; Frequency band to use (Hz). If you 
%                        choose a different band, consider changing 
%                        cfg.kernel accordingly.
%       cfg.smooth = 1; If 1 apply smoothing, if 0 don't. Smoothing is
%                       recommended for better performance.
%       cfg.kernel = []; Kernel to use for smoothing. The default kernel is 
%                        gausskernel(60,20), and is used when cfg.kernel is
%                        empty. cfg.kernel should be [n x 1] double.
%       cfg.verbose = 1; Tell me what you are doing.   
%
% (proposed mundane name for OldWizard: htSWR for "hilbert transform SWR")
% aacarey Nov 2015.
%
%  see also TSDtoIV2

%% Parse cfg parameters
cfg_def.rippleband = [140 250]; % in Hz
cfg_def.smooth = 1; % do you want to smooth the detector or not
cfg_def.kernel =[]; % which kernel (if empty, goes to default)
cfg_def.verbose = 1; % talk to me or not

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose
    tic
    disp([mfun,': looking for sharp wave-ripple events...'])
end

% filter in the ripple band

cfg_temp = [];
cfg_temp.type = 'fdesign';
cfg_temp.f = cfg.rippleband;
CSCf = FilterLFP(cfg_temp,CSC);

% ask hilbert what he thinks
score = abs(hilbert(CSCf.data)).^2;

% apply smoothing, if desired

if cfg.smooth
    if isempty(cfg.kernel)
        kernel = gausskernel(60,20);
    else
        kernel = cfg.kernel;
    end   
   score = conv(score,kernel,'same'); 
end

% make output
SWR = tsd(CSC.tvec,score);

SWR.cfg.history.mfun = cat(1,SWR.cfg.history.mfun,mfun);
SWR.cfg.history.cfg = cat(1,SWR.cfg.history.cfg,{cfg});

if cfg.verbose
    toc
end

end

