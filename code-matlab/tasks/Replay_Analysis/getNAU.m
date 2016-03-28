function NAU = getNAU(cfg_in,S);
%% GETNAU Get Number of Active Units
%   NAU = getNAU(cfg_in,S) returns a tsd with number of independantly active units
%
%   INPUT:
%       cfg_in - input cfg parameters
%       S - input ts with spiketrains
%
%   OUTPUT:
%       NAU - tsd with number of active units
%
%   CFG OPTIONS:
%       cfg.tvec = []; % Internal time vector
%       cfg.sigma = 0.05;
%       cfg.dt = 0.001; % internal timebin size in seconds
%       cfg.convFlag = 0;
%       cfg.removenan = 1;
%
% youkitan 2015-01-01

%% Parse cfg parameters
cfg_def.tvec = [];
cfg_def.sigma = 0.05;
cfg_def.dt = 0.001; % internal timebin size in seconds
cfg_def.convFlag = 0;
cfg_def.removenan = 1;

cfg = ProcessConfig2(cfg_def,cfg_in); % this takes fields from cfg_in and puts them into cfg
mfun = mfilename;

%% Get spikes and tvec

% concatenate all spikes
spk = [];
for iC = 1:length(S.t)
    spk = cat(1,spk,S.t{iC});
end

% construct internal tvec
if ~isempty(cfg.tvec)
    if ~cfg.convFlag
        tveci_edges = cfg.tvec;
    else
        tveci_edges = cfg.tvec(1):cfg.dt:cfg.tvec(end);
    end
    
else
    sprintf('warning: cfg.tvec is not defined or empty. Making new internal tvec!');
    tveci_edges = min(spk)-cfg.dt:cfg.dt:max(spk)+cfg.dt;
end

tveci_centers = tveci_edges(1:end-1)+cfg.dt/2;

%% Make NAU tsd
NAU = tsd; NAU.label{1} = 'NAU';

% Count spikes per bin, for each cell
spk_binned = zeros(length(S.t),length(tveci_centers));
for iC = 1:length(S.t)
    spk_hist = histc(S.t{iC},tveci_edges)';
    spk_hist(end-1) = spk_hist(end-1)+spk_hist(end); % combine last bin and edge
    spk_hist = spk_hist(1:end-1);
    spk_binned(iC,:) = spk_hist;
end

% count number of cells firing in each bin 
if cfg.convFlag %convolution of counts
    % construct gaussian kernel
    gauss_window = 1./cfg.dt; 
    gauss_SD = cfg.sigma./cfg.dt; % 0.02 seconds (20ms) SD
    gk = gausskernel(gauss_window,gauss_SD); 
    gk = gk./cfg.dt; % normalize by binsize
    
    spk_conv = zeros(length(S.t),length(tveci_centers));
    
    %for each cell, convolve the spike train
    for iC = 1:length(S.t)
        spk_conv(iC,:) = conv(spk_binned(iC,:),gk,'same');
    end
    
    spk_out = sum(spk_conv);
    
    % interpolate back on original timebase
    NAU.data = interp1(tveci_centers,spk_out,cfg.tvec,'linear')';
    NAU.tvec = cfg.tvec;

else % pure counts
    spk_out = sum(spk_binned ~= 0);
    NAU.data = spk_out;
    NAU.tvec = tveci_centers;
    
    NAU2.data = interp1(tveci_centers,spk_out,cfg.tvec,'linear')';
    NAU2.tvec = cfg.tvec;
end

% remove nans if desired
if cfg.removenan
    keep_idx = ~isnan(NAU.data);
    NAU.data = NAU.data(keep_idx); NAU.tvec = NAU.tvec(keep_idx);
end

% housekeeping
cfg = rmfield(cfg,'tvec'); % is in MUA itself anyway
NAU.cfg.history.mfun = cat(1,NAU.cfg.history.mfun,mfun);
NAU.cfg.history.cfg = cat(1,NAU.cfg.history.cfg,{cfg});