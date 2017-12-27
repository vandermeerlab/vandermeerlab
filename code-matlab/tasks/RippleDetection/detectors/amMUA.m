function [MUA,raw,noise] = amMUA(cfg_in,S,tvec)
%AMMUA Regions of concentrated multi-unit spiking activiy
%   [MUA,raw,noise] = amMUA(cfg_in,S,tvec) uses adaptive thresholding to identify regions 
%   containing greater multi-unit activity than neighbouring regions. 
%
%   The defaults are recommended. Change one config at a time and plot MUA,
%   noise, and raw curves along with S to decide on the best parameters. 
%
%   Pretty slow.
%
%   INPUTS
%   S    - output from LoadSpikes.
%   tvec - time vector from LoadCSC 
%
%   CONFIGS (with defaults):
%
%   cfg.spkcap = 2; Factor corresponding to how much weight each cell can contribute.
%                   If a cell spikes more than this number of times in
%                   rapid succession, it only gets counted as this many
%                   spikes. The region containing the capped spikes spans
%                   the standard deviation of the kernel. Increasing this
%                   number increases the allowed "loudness" of a cell. 
%                   If spkcap is 0, the curve is flat.
%
%   cfg.noisefloor = 4; Magic constant for how much noise we try to remove.
%                    Measured in # of cells. How many noisy cells will we
%                    tolerate? Increasing this number narrows the detection
%                    width for individual events.
%   cfg.verbose = 1; talk to me
%
% OUTPUTS
% MUA - tsd with fields
%       .data: curve identifying regions of multi-unit 
%          activity 
%       .tvec: pointless copy of csc.tvec
%
% raw - tsd
%       .data: raw muascore before adaptive thresholding
%       .tvec: time vector
%
% noise - tsd
%       .data: curve used for adaptive thresholding; raw - noise = MUA
%       .tvec: time vector
%
% Elyot Grant, Feb 2015 (adaptive thresholding upgrade)
% ACarey, Jan 2015 

%% Parse cfg parameters

cfg_def.spkcap = 2;
cfg_def.noisefloor = 4; 
cfg_def.kernelstd = 40; 
cfg_def.verbose = 1;
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

%%

if isempty(S.t) 
    error('No spikes: cannot perform detection.')
end

if cfg.verbose
    tic
    disp([mfun,': Looking for multi-unit activity...']);
    disp(' ');
end

%% muascore: regions of higher multi-unit activity
%tic

%skpcap = 2; %Factor corresponding to how much weight each cell can contribute.
                    %If a cell spikes more than this number of times in
                    %rapid succession, it only gets counted as this many
                    %spikes.
 kernelstd = cfg.kernelstd;                   
cellactivitycap = cfg.spkcap./kernelstd./sqrt(2*pi);

muascore = zeros(size(tvec));
%muascoreneg = zeros(size(tvec));
%muascoretot = zeros(size(tvec));
kernel = gausskernel(5*kernelstd,kernelstd);
%kernel2 = -5*gausskernel(3000,300);
for iS = 1:length(S.t)
    spk_tvec = zeros(size(tvec));
    spiketrain = S.t{iS};
    for iSpike = 1:length(spiketrain)
        spk_here = nearest_idx3(spiketrain(iSpike),tvec);
        spk_tvec(spk_here) = spk_tvec(spk_here)+1;
    end
    spk_conv = conv(spk_tvec,kernel,'same');
    %spk_conv2 = conv(spk_tvec,kernel2,'same');
    spk_conv = min(spk_conv,cellactivitycap);
    muascore = muascore + spk_conv;
    %muascoreneg = muascoreneg + 0.1*spk_conv2;
    
    %spk_tot = spk_conv + spk_conv2;
    %spk_tot = max(0,spk_tot);
    %spk_tot = min(spk_tot,cap);
    %muascoretot = muascoretot + spk_tot;
end

% Rescale
themean = mean(muascore);

%cfg.noisefloor = 4; %Magic constant for how much noise we try to remove.
                %measured in # of cells.
weightednoisefloor = cellactivitycap.*cfg.noisefloor;
detectionthreshold = 1; %How many cells needed for detection (above noise floor)
% ^^^ this is stupid, should not be hardcoded...not to mention there's a
% variable name assigned to the number 1.

muacapped = min(weightednoisefloor,muascore);
kernel3 = -gausskernel(3000,250);
muascoreneg2 = conv(muacapped,kernel3,'same');

 muascore = muascore ./ themean;
%muascoreneg = 3.5*muascoreneg ./ themean;
 muascoreneg2 = muascoreneg2 ./ themean;
 %muascoretot = muascoretot ./ themean;
 %muascoretot2 = max(0,muascore + muascoreneg);
 muascoretot3 = max(0,muascore + muascoreneg2 - detectionthreshold.*cellactivitycap ./ themean);

%MUA2 = tsd(tvec,muascoreneg);

%MUA3 = tsd(tvec,muascoretot);
%MUA = tsd(tvec,muascoretot2); % old output

%MUA5 = tsd(tvec,muascoretot2.*muascore/mean(muascore));

%% Generate output

MUA = tsd(tvec,muascoretot3'); % raw-noise
MUA = History(MUA,mfun,cfg);

raw = tsd(tvec,muascore'); % raw
raw = History(raw,mfun,cfg);

noise = tsd(tvec,muascoreneg2'); % noise
noise = History(noise,mfun,cfg);

if cfg.verbose
toc
end
end

