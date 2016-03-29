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
%
%       cfg.weightby = 'amplitude'; Do you want the amplitude envelope or
%                       the power envelope?
%            'amplitude' - abs(hilbert(filtered_signal))
%            'power'     - abs(hilbert(filtered_signal)).^2
%
%       cfg.smooth = 1; If 1 apply smoothing, if 0 don't. Smoothing is
%                       recommended for better performance.
%
%       cfg.kernel = 'wizard'; Kernel to use for smoothing. See cfg.kernel
%                        options in ConvTSD.
%
%       cfg.verbose = 1; If 1, print informative text to the command
%                        window, if 0, don't.
%
% (proposed mundane name for OldWizard: htSWR for "hilbert transform SWR")
% aacarey Nov 2015.
%
%  see also TSDtoIV2

%% Parse cfg parameters
cfg_def.rippleband = [140 250]; % in Hz
cfg_def.weightby = 'amplitude';
cfg_def.smooth = 1; % do you want to smooth the detector or not
cfg_def.kernel = 'wizard'; % which kernel (if empty, goes to default)
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
cfg_temp.verbose = 0;
CSCf = FilterLFP(cfg_temp,CSC);

% ask hilbert what he thinks

switch cfg.weightby
    case 'amplitude'
        score = abs(hilbert(CSCf.data));
    case 'power'
        score = abs(hilbert(CSCf.data)).^2;
    otherwise
        error('Unrecognized option specified in cfg.weightby')
end

% convert to tsd struct
SWR = tsd(CSC.tvec,score);

% apply smoothing, if desired
if cfg.smooth
    cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.kernel = cfg.kernel; cfg_temp.where = 'all';
    SWR = ConvTSD(cfg_temp,SWR);
end

SWR = History(SWR,mfun,cfg);

if cfg.verbose
    toc
end

end

