%% Parameters
dsRate = 20; % Downsampling rate if you want to downsample the signal
%This dsRate will also be applied to all signals during the analysis
%pipeline

% Filter Parameters
lpCut = 20; % Cut-off frequency for filter
filtOrder = 8; % Order of the filter

% Baseline Parameters
basePrc = 5; % Percentile value from 1 - 100 to use when finding baseline points
%Note: Lower percentiles are used because the mean of signal is not true
%baseline
winSize = 10; % Window size for baselining in seconds
winOv = 0; %Window overlap size in seconds
interpType = 'linear'; % 'linear' 'spline' 
fitType = 'interp'; % Fit method 'interp' , 'exp' , 'line'

% Demodulation Parameters
%When demodulating signals, the filter creates edge artifacts. We record
%for a few seconds longer, so we can remove x seconds from the beginning
%and end
%Adjust the variable to "0" if it's a normal photometry recording
sigEdge = 15; %Time in seconds of data to be removed from beginning and end of signal
modFreq = [319 217];

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

%% Data preparation
rawFs = FP_data.acq.Fs;
Fs = rawFs/dsRate;
rawFP = FP_data.acq.FP{1}';
rawTime = FP_data.acq.time';
n_timeVec = size(rawTime, 2);
n_timeVec = length(1:dsRate:n_timeVec);
timeVec = [1:n_timeVec]/Fs;

%% Plotting raw data
sessionTitle = 'M21-061R_';
time_ranges = [10, 50, 100];

for t_i = 1:length(time_ranges)
    t_range = 1:rawFs*time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(rawTime(t_range), rawFP(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fluorescence (dF/F)'); xlabel('Time (s)');
end

%% Filtering and downsampling
nbFP = filterFP(rawFP,rawFs,lpCut,filtOrder,'lowpass');
nbFP = downsample(nbFP, dsRate);

for t_i = 1:length(time_ranges)
    t_range = 1:Fs*time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(timeVec(t_range), nbFP(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fluorescence (dF/F)'); xlabel('Time (s)');
end

%% Prepare baseline
%Ensure the FP vector is a column vector instead of row vector --> Faster
%computation
if size(nbFP,1) == 1
    nbFP = nbFP';
end

winSize = winSize * Fs; %Convert window size from seconds to samples
winOv = winOv * Fs; %Convert overlap window from seconds to samples
Ls = length(nbFP); %Get length of photometry trace
L = 1:Ls; L = L'; %Create a column vector from 1 to total data points in trace
nPts = floor(Ls/(winSize-winOv)); %Determine number of baseline points to be found

%X is a vector of positional points in the data that we will be gathering baseline points
%Y is a vector of zeros that will contain calculated baseline points
X = L(1:ceil(Ls/nPts):end); Y = zeros(nPts,1);

%Determine the step size of the window:
%If the overlap is 0 or empty then it will use the window size as the step
%size. If the overlap is greater than 0 the step size will be the window
%size subtracted by the overlap size
if winOv == 0 || isempty(winOv)
    winStep = winSize;
else
    winStep = winSize - winOv;
end

%The following for loop goes through the photometry vector and finds
%baseline values of the windowed photometry trace according to a certain
%percentile
for n = 0:nPts-1
    I1 = (n*winStep)+1;
    I2 = I1 + winSize;
    if I2>Ls
        I2 = Ls;
    end
    Y(n+1) = prctile(nbFP(I1:I2),basePrc);
end

interpFit = interp1(X,Y,L,interpType,'extrap'); %Create an interpolated line using the previously calculated baseline points
baseline = interpFit; %Use interpolated line

%% Plotting baseline
for t_i = 1:length(time_ranges)
    t_range = 1:Fs*time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(timeVec(t_range), baseline(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fluorescence (dF/F)'); xlabel('Time (s)');
end

%% Detrend from baseline
dF = (nbFP-baseline)./baseline;
dF = dF*100;

%% Plotting dF/F (normalized fluorescence level)
for t_i = 1:length(time_ranges)
    t_range = 1:Fs*time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(timeVec(t_range), dF(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fluorescence (dF/F)'); xlabel('Time (s)');
end