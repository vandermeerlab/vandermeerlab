function [SWR,swr1,swr2] = amSWR(cfg_in,SWRfreqs,csc)
%AMSWR Detect putative sharp wave-ripple events using the discrete 
% Fourier transform.
%
% function SWR = amSWR(cfg_in,SWRfreqs,csc)
%
%   INPUTS
%   SWRfreqs - Fourier coefficient struct for SWRs (output from SWRfreak).
%   csc - output from LoadCSC, comprising the data from a continuously sampled 
%         channel (csc). 
%
%   CONFIGS (with defaults):
%   cfg.fs = 2000; in Hz, the sampling frequency
%
%   OUTPUTS
%   SWR - [nSamples x 1 double] timestamped data struct 
%              containing unitless values that describe how closely csc.data 
%              resembled an SWR. 
%  
% Elyot Grant, Jan 2015 (detection)
% ACarey, Jan 2015 (function design)

%% Parse cfg parameters

cfg_def.fs = 2000;
cfg_def.verbose = 1;
cfg = ProcessConfig2(cfg_def,cfg_in);

%% Tell the user you're doing something


if cfg.verbose
    tic
    cprintf(-[0 0 1],'amSWR: Identifying potential sharp wave-ripple events...');
    disp(' ');
end

 %% Set up 
 
 %freqs = [0,0,0,0,0,0,0,0,0,0,0,0,0.000170521028425372,0.00170778248757081,0.00331395105024794,0.00479391795042840,0.00607136317217047,0.00678590934283586,0.00697574495989578,0.00650036253423614,0.00569838612070232,0.00466283089230949,0.00366655000685789,0.00272766671935131,0.00195939311675901,0.00137725227252137,0.00102976121026613,0.000775147516298750,0.000613151554386769,0.000445041806482364,0.000371607082692109,0.000310870018623868,0.000321967274810346,0.000312497164011530,0.000308075500895654,0.000273171891686851,0.000204040282703605,0.000109334819915455,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
 
    function swrscore = SWRhelper(freqs,win,csc,fs)
 
        sampwin = win*fs; % the window size in nSamples = win * sampling frequency
 
        hiPassCutoff = 100; %We want to delete all frequencies below 100 Hz
        fourierCoeffCutoff = round(hiPassCutoff*win);
 
        swrscore = zeros(length(csc.data),1);
 
        % FFT
        %for iSWR = sampwin:500000-sampwin % for testing
        for iSWR = sampwin:length(csc.data)-sampwin
        %for iSWR = 3868190:3868191
    
            % Grab a windowed candidate and smooth the ends to reduce noise.
            next = windowedFFT(csc.data,fs,sampwin,iSWR);
            next(1:fourierCoeffCutoff) = 0; % delete all frequencies below 100 Hz (this is only for norming)
            %nextNormed = next / sum(next);

            score = sum(next.*freqs);
            %score2 = sum(nextNormed.*freqs);
            swrscore(iSWR) = score;

        end
        
        % get rid of imaginary numbers (very rare, but happens)
        swrscore = max(0,swrscore); 
    end

%% get score vectors

 swrscore1 = SWRhelper(SWRfreqs.freqs1,SWRfreqs.win1,csc,cfg.fs);
 % rescale because values are insanely low
 swrscore1 = rescmean(swrscore1,1);
 
 disp('*yawn* ...')
 
 swrscore2 = SWRhelper(SWRfreqs.freqs2,SWRfreqs.win2,csc,cfg.fs);
 % rescale
 swrscore2 = rescmean(swrscore2,1);

 % combine the scores
 geometricmean = sqrt(swrscore1.*swrscore2);
 
 %% return
 SWR = tsd(csc.tvec,geometricmean);
 SWR.win1 = SWRfreqs.win1;
 SWR.win2 = SWRfreqs.win2;
 swr1 = tsd(csc.tvec,swrscore1);
 swr2 = tsd(csc.tvec,swrscore2);
 
 if cfg.verbose
 toc
 end
end