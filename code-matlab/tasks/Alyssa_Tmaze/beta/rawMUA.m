function [MUA] = rawMUA(cfg_in,S,tvec)
%RAWMUA Raw multi-unit activity from tt files
% [MUA] = RAWMUA(cfg_in,S,tvec) returns a TSD containing raw multi-unit
%   activity. The user can choose to smooth the output (this is default,
%   recommended).
%   
%    INPUTS
%       S: ts datatype containing spike times, defaults assume S is from
%          .tt files from dCA1 (lots of spikes!) so you may need to use a
%          custom kernel or see getMUA for .t file spikes.
%
%       tvec: the time vector that spikes will be binned into, ex: CSC.tvec
%
%    OUTPUTS
%        MUA: TSD containing binned, smoothed raw multi-unit activity (raw
%             as in, not adaptively rescaled)
% 
%    CONFIG OPTIONS
%        cfg.smooth = 1; If 1, smooth the output, if 0 don't.
%
%        cfg.kernel = 'wizard'; Applies smoothing with selected kernel if 
%                          cfg.smooth = 1. See cfg.kernel in ConvTSD for
%                          options. ('wizard' and 'fang' seem good)
%
% see also getMUA, MakeTTFiles
%
% aacarey Aug 2015, edit Jan 2015, Feb 2015

%% Set config defaults and parse parameters

cfg_def.verbose = 1;
cfg_def.smooth = 1;
cfg_def.kernel = 'wizard';

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% Put all spikes into a single vector and sort them
allSpikes = [];
for iTrain = 1:length(S.t)
   allSpikes = cat(1,allSpikes,S.t{iTrain}); 
end

allSpikes = sort(allSpikes);

% bin the spikes into the time vector using nearest_idx3...this makes sense
% only if the tvec and S have been restricted to the same times. If S is
% unrestricted, some spikes will be binned incorrectly.
spikebin = zeros(size(tvec)); % spikebin in raw multiunit activity
for iSpike = 1:length(allSpikes)
    spike_here = nearest_idx3(allSpikes(iSpike),tvec);
    spikebin(spike_here) = spikebin(spike_here)+1;
end

% spikebin has collected all the spikes. Some bins will have multiple
% spikes, and some will have zero.

% make output
MUA = tsd(tvec,spikebin');

% Apply smoothing
if cfg.smooth
    cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.where = 'all'; cfg_temp.kernel = cfg.kernel;
    MUA = ConvTSD(cfg_temp,MUA);
end

% make output
MUA = History(MUA,mfun,cfg);

end

