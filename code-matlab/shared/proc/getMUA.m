function MUA = getMUA(cfg_in,S)
% function MUA = getMUA(cfg,S)
%
% INPUT:
% ts with spike trains
%
% OUTPUT:
% tsd with multi-uinit activity
%
% cfg options with defaults:
%
% cfg.sigma = 0.05; % standard deviation of gaussian for convolution
% cfg.tvec = []; % MUST BE SPECIFIED: timebase for output tsd
% cfg.dt = 0.001; % internal timebase
% cfg.removenan = 1;
%
% MvdM 2014-06-25

cfg_def.sigma = 0.05;
cfg_def.tvec = [];
cfg_def.dt = 0.001; % internal timebase
cfg_def.removenan = 1;

cfg = ProcessConfig2(cfg_def,cfg_in); % this takes fields from cfg_in and puts them into cfg
mfun = mfilename;

if ~isfield(cfg,'tvec')
   error('cfg.tvec must be defined.'); 
end

MUA = tsd; MUA.label{1} = 'MUA';

% construct internal tvec
tveci_edges = cfg.tvec(1):cfg.dt:cfg.tvec(end);
tveci_centers = tveci_edges(1:end-1)+cfg.dt/2;

% concatenate spikes into one list
spk = [];
for iC = 1:length(S.t)
   spk = cat(1,spk,S.t{iC}); 
end

% bin list
spk_binned = histc(spk,tveci_edges);
spk_binned(end-1) = spk_binned(end-1)+spk_binned(end); % combine last bin and edge
spk_binned = spk_binned(1:end-1);

% construct gaussian kernel
gauss_window = 1./cfg.dt; 
gauss_SD = cfg.sigma./cfg.dt; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg.dt; % normalize by binsize

% convolve
spk_conv = conv2(spk_binned,gk,'same');

% interpolate back on original timebase
MUA.data = interp1(tveci_centers,spk_conv,cfg.tvec,'linear')';
MUA.tvec = cfg.tvec;

% remove nans if desired
if cfg.removenan
    keep_idx = ~isnan(MUA.data);
    MUA.data = MUA.data(keep_idx); MUA.tvec = MUA.tvec(keep_idx);
end

% housekeeping
cfg = rmfield(cfg,'tvec'); % is in MUA itself anyway
MUA.cfg.history.mfun = cat(1,MUA.cfg.history.mfun,mfun);
MUA.cfg.history.cfg = cat(1,MUA.cfg.history.cfg,{cfg});