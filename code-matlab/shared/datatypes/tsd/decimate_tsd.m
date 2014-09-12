function tsd_in = decimate_tsd(cfg_in,tsd_in)
% function tsd_in = decimate_tsd(cfg,tsd_in)
%
%
% DOES NOT YET WORK FOR MULTIDIMENSIONAL DATA

cfg = [];
cfg.decimateFactor = 4;
ProcessConfig;

mfun = mfilename;

tsd_in.data = decimate(tsd_in.data,cfg.decimateFactor);
tsd_in.tvec = downsample(tsd_in.tvec,cfg.decimateFactor);

tsd_in.cfg.hdr{1}.SamplingFrequency = tsd_in.cfg.hdr{1}.SamplingFrequency ./ cfg.decimateFactor;

% housekeeping
tsd_in.cfg.history.mfun = cat(1,tsd_in.cfg.history.mfun,mfun);
tsd_in.cfg.history.cfg = cat(1,tsd_in.cfg.history.cfg,{cfg});