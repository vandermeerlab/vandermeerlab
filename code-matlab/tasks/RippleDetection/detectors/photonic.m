function tsd_out = photonic(cfg_in,tsd_in)
%PHOTONIC Detect transients in the local field potential
%
%  TSD_out = PHOTONIC(cfg,TSD)
%
%    Transient detection works best with little background noise. Since the 
%  HC has a lot of high-amplitude, low frequency background "noise", you 
%  should highpass or bandpass filter the signal before sending it in (see 
%  FilterLFP and fftFilterLFP).
%
%    INPUTS
%          cfg: config struct with fields controlling function behavior
%          TSD: tsd struct
%
%    OUTPUTS
%      TSD_out: tsd struct with .data field containing detection results
%
%    CONFIG OPTIONS 
%       cfg.weightby = 'power'; How to obtain the signal envelope
%          'amplitude'    - amplitude envelope (absolute value of the signal)
%          'power'        - power envelope (square the signal)
%
%       cfg.useHilbert = 1; If 1, uses Hilbert transform on input, if 0
%                     doesn't
% 
%       cfg.kernel = 'wizard'; Which smoothing kernel to use. See 
%                     cfg.kernel in ConvTSD.
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
cfg_def.weightby = 'amplitude';
cfg_def.useHilbert = 0;
cfg_def.kernel = 'wizard';
cfg_def.smooth = 1;
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
        if cfg.useHilbert
            envelope = abs(hilbert(tsd_in.data));
        else
            envelope = abs(tsd_in.data);
        end
    case 'power' % power envelope
        if cfg.useHilbert
            envelope = abs(hilbert(tsd_in.data)).^2;
        else
            envelope = tsd_in.data.^2;
        end
    otherwise
        error('Unrecognized config option specified in cfg.weightby')
end

tsd_out = tsd(tsd_in.tvec,envelope);

if cfg.smooth
    cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.where = 'all'; cfg_temp.kernel = cfg.kernel;
    tsd_out = ConvTSD(cfg_temp,tsd_out);
end

% write config history
tsd_out = History(tsd_out,mfun,cfg);

end

