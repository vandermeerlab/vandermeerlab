function fftedData = windowedFFT(cfg,data,sampwin,idx)
%WINDOWEDFFT helper function for SWRfreak and amSWR

%   Elyot Grant, Jan 2015

fs = cfg.fs;
hiPassCutoff = cfg.hiPassCutoff;

%hiPassCutoff = 100; %We want to delete all frequencies below 100 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSmoothingSamples = round(1*fs/hiPassCutoff); % Number of samples in the wavelength (frequency) specified by hiPassCutoff
newSampwin = sampwin - numSmoothingSamples;
smoothingWin1 = 0.5 - 0.5*cos(pi/numSmoothingSamples*(0:numSmoothingSamples-1));
smoothingWin2 = ones(1,newSampwin);
smoothingWin3 = 0.5 - 0.5*cos(pi/numSmoothingSamples*(numSmoothingSamples-1:-1:0));
smoothingWin = [smoothingWin1 smoothingWin2 smoothingWin3];

windowedData = data(idx-floor(newSampwin/2)-numSmoothingSamples:idx+ceil(newSampwin/2)+numSmoothingSamples-1);

windowedData = windowedData .* smoothingWin; % Smooth the window.

% Overlap the window with itself
windowedData(1:numSmoothingSamples) = windowedData(1:numSmoothingSamples) + windowedData(newSampwin+numSmoothingSamples+1:newSampwin+2*numSmoothingSamples);
windowedData = windowedData(1:numSmoothingSamples+newSampwin);

fftedData = abs(fft(windowedData)); % abs converts complex numbers to their magnitudes (because we don`t care about the phase)
fftedData = fftedData(1:round(length(fftedData)/2));

%% "weight by power", but not quite...should remove this option and stop using the terms "amplitude-weighted" and "power-weighted" 
% divide by 1/x...aka multiply by x...
if strcmp(cfg.weightby,'power')
    fftedData = fftedData.*(1:length(fftedData)); % weight by power
end
