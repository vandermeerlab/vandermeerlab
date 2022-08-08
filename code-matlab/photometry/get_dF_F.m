%% Parameters
params.dsRate = 20; % Downsampling rate if you want to downsample the signal
%This dsRate will also be applied to all signals during the analysis
%pipeline

% Filter Parameters
params.FP.lpCut = 20; % Cut-off frequency for filter
params.FP.filtOrder = 8; % Order of the filter

% Baseline Parameters
params.FP.basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline
params.FP.winSize = 10; % Window size for baselining in seconds
params.FP.winOv = 0; %Window overlap size in seconds
params.FP.interpType = 'linear'; % 'linear' 'spline' 
params.FP.fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'

% Demodulation Parameters
%When demodulating signals, the filter creates edge artifacts. We record
%for a few seconds longer, so we can remove x seconds from the beginning
%and end
%Adjust the variable to "0" if it's a normal photometry recording
params.FP.sigEdge = 15; %Time in seconds of data to be removed from beginning and end of signal
params.FP.modFreq = [319 217];

%%
addpath(genpath('/Users/mac/Projects/nsb2021/code-matlab/shared'));
addpath(genpath('/Users/mac/Projects/nsb2021/photometry'));

%% Load Neuralynx CSC photometry data
cfg.fc = {'CSC33.ncs'};
csc_photo = LoadCSC(cfg);

FP_data = [];
FP_data.acq.Fs = csc_photo.cfg.hdr{1}.SamplingFrequency;
FP_data.acq.time = csc_photo.tvec - csc_photo.tvec(1);
FP_data.acq.FP{1} = csc_photo.data';

%% Process data recorded in CW mode
FP_data = processFP(params, FP_data);

FP = tsd(FP_data.final.time', FP_data.final.FP{1}', 'FP');

%% Load events
LoadExpKeys();
cfg_evt = [];
cfg_evt.eventList = ExpKeys.eventList;
cfg_evt.eventLabel = ExpKeys.eventLabel;
evt = LoadEvents(cfg_evt);

%% PETH
left_end = evt.t{3};
right_end = evt.t{4};

cfg_peth.dt = 0.01; % time step, specify this for 'interp' mode
cfg_peth.mode = 'interp'; % 'raw' or 'interp'; if 'interp', need to specify cfg_def.dt
left_end_peth = TSDpeth(cfg_peth, FP, left_end);
right_end_peth = TSDpeth(cfg_peth, FP, right_end);