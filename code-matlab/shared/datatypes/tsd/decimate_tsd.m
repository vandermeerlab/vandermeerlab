function tsd_in = decimate_tsd(cfg_in,tsd_in)
% function tsd_in = decimate_tsd(cfg,tsd_in)
%
%
% DOES NOT YET WORK FOR MULTIDIMENSIONAL DATA

cfg = [];
cfg.decimateFactor = 4;
cfg = ProcessConfig(cfg, cfg_in);

mfun = mfilename;

if ~CheckTSD(tsd_in)
   error('Input data is not correctly formed.');
end

for iData = size(tsd_in.data,1):-1:1
    
    data_temp(iData,:) = decimate(tsd_in.data(iData,:),cfg.decimateFactor);
    tsd_in.cfg.hdr{iData}.SamplingFrequency = tsd_in.cfg.hdr{iData}.SamplingFrequency ./ cfg.decimateFactor;
    
end

tsd_in.data = data_temp;
tsd_in.tvec = downsample(tsd_in.tvec,cfg.decimateFactor);

% housekeeping
tsd_in.cfg.history.mfun = cat(1,tsd_in.cfg.history.mfun,mfun);
tsd_in.cfg.history.cfg = cat(1,tsd_in.cfg.history.cfg,{cfg});