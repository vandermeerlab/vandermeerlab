%% WORKFLOW_process_photometry_data
% basic script to plot raw photometry data acquired with Neuralynx,
% compute deltaF/F, and plot the result

%% set path 
addpath('C:\Users\mvdm\Documents\GitHub\nsb2022\Photometry'); % in addition to shared code

%% parameters
params.dsRate = 20; % Downsampling rate (in Hz) if you want to downsample the signal

% Filter Parameters
params.FP.lpCut = 20; % Cut-off frequency for filter
params.FP.filtOrder = 8; % Order of the filter

% Baseline Parameters (used to compute deltaF/F)
params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline
params.FP.winSize = 10; % Window size for baselining in seconds
params.FP.winOv = 0; %Window overlap size in seconds
params.FP.interpType = 'linear'; % 'linear' 'spline' 
params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'

% Demodulation Parameters (when using isobestic)
%When demodulating signals, the filter creates edge artifacts. We record
%for a few seconds longer, so we can remove x seconds from the beginning
%and end
%Adjust the variable to "0" if it's a normal photometry recording
params.FP.sigEdge = 0; %Time in seconds of data to be removed from beginning and end of signal
params.FP.modFreq = [319 217]; % modulation frequencies for signal and isobestic light sources

%% Load Neuralynx CSC photometry data
cfg = [];
cfg.fc = {'CSC33.ncs'};
csc_photo = LoadCSC(cfg);

subplot(211)
plot(csc_photo);
set(gca, 'FontSize', 18, 'LineWidth', 1, 'TickDir', 'out');
title('raw data'); box off; axis tight;

%% Process data recorded in CW mode (default)
FP_data = [];
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency;
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1);
FP_data.acq.FP{1} = csc_photo.data';

FP_data = processFP(params, FP_data); %% <-- this is the key processing step
FP = tsd(FP_data.final.time', FP_data.final.FP{1}', 'FP');

%% plot
subplot(212)
plot(FP);
set(gca, 'FontSize', 18, 'LineWidth', 1, 'TickDir', 'out');
title('preprocessed data'); box off;
ylabel('\DeltaF/F (%%)'); xlabel('time (s)'); axis tight;