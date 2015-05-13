function SWRfreqs = SWRfreak(cfg_in,SWRtimes,csc)
%SWRFREAK Fourier coefficients of summed sharp wave-ripples.
%
%   SWRfreqs = swrfreak(SWRtimes,csc,fs)
%
%   SWRfreqs = swrfreak(SWRtimes,csc,fs,NOISEtimes)
%
%   NOTE: you can change trapwin at this step, but you must use the same
%   one in detectSWR
%
%   cfg.win1 = 0.08; in ms, the "noise reduction" window
%   cfg.win2 = 0.04; in ms, the "precision" window
%   cfg.fs = 2000; in Hz, the sampling frequency
%
%   OUTPUT
%   
%   SWRfreqs - frequency data for SWRs, in terms of Fourier coefficients.
% 
%   For additional information, open SWRfreak and read the ABOUT section
%
% Elyot Grant, Jan 2015 
% ACarey, Jan 2015

%% ABOUT SWRFREAK

% Why are there two windows, and why are 40 ms and 80 ms the defaults?

% The noise reduction and precision windows were chosen based on how well
% the plotted scores (from amSWR) characterized the data in the csc and the 
% spiketrains.
% 120, 100, 80, 60, and 40 ms windows were used for FFTing the data. 120
% and 100 merged nearby events and were rejected, despite distinguishing
% between false positives better. 80 ms still merged nearby events, but to
% a lesser degree than 100 and 120. It also had pretty good noise
% reduction. 40 ms had the best separation for nearby events, but was the
% worst for peaking at false positives. The geometric mean of the two of
% them did a really good job, so it was decided that two sets of SWRfreqs
% were better than one. 

%% Parse cfg parameters

cfg_def.win1 = 0.08;
cfg_def.win2 = 0.04;
cfg_def.fs = 2000; 

cfg = ProcessConfig2(cfg_def,cfg_in);

%% check if csc is the same as the one used in ducktrap

if ~strcmp(SWRtimes.label,csc.label)
    warning('CSC is different from the one used for manual SWR identification')
end

%% internal function to do the thing

    function freqs = freakHelper(fs,SWRtimes,csc,timewin)
        sampwin = timewin*fs; % the window size in nSamples = timewindow * sampling frequency 
        %tstart = SWRtimes.tcent - timewin/2;
        %tend = SWRtimes.tcent + timewin/2;
        midIndices = nearest_idx3(SWRtimes.tcent,csc.tvec);
        SWRsum = 0;
        for iSWR = 1:length(SWRtimes.tcent) 
            nextFFT = windowedFFT(csc.data,fs,sampwin,midIndices(iSWR));
            SWRsum = SWRsum + nextFFT;
        end

        midIndices = nearest_idx3(SWRtimes.tcent+2,csc.tvec); % Add 2 seconds to each clicked time -> random time (this could error if a clicked time was closer than 2 s away from the end of recording?)
        SWRsumNoise = 0;
        for iSWR = 1:length(SWRtimes.tcent)
            nextFFT = windowedFFT(csc.data,fs,sampwin,midIndices(iSWR));
            SWRsumNoise = SWRsumNoise + nextFFT;
        end

        SWRsum = conv(SWRsum,[0.1,0.2,0.4,0.2,0.1]); % narrow smoothing kernel
        SWRsum = SWRsum(3:length(SWRsum)-2);
        SWRsum = SWRsum./sum(SWRsum); %Normalize
        %plot(SWRsum);

        SWRnoise = conv(SWRsumNoise,[0.1,0.2,0.4,0.2,0.1]);
        SWRnoise = SWRnoise(3:length(SWRnoise)-2);
        SWRnoise = SWRnoise./sum(SWRnoise);
        %figure;plot(SWRnoise);

        freqs = (SWRsum - SWRnoise);
        freqs = conv(freqs,[0.1,0.2,0.4,0.2,0.1]);
        freqs = freqs(3:length(freqs)-2);
        %freqs = max(freqs,zeros(size(freqs))); %%%%%%%%%%%%%%%%%%%%%%
        %figure;plot(SWRfreqs);
    end

%% generate the output

freqs1 = freakHelper(cfg.fs,SWRtimes,csc,cfg.win1);
freqs2 = freakHelper(cfg.fs,SWRtimes,csc,cfg.win2);

%% return the output
SWRfreqs = struct('freqs1',freqs1,'freqs2',freqs2,'win1',cfg.win1,'win2',cfg.win2);
SWRfreqs.label = csc.label;

end

