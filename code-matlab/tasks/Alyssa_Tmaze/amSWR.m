function [SWR,swr1,swr2] = amSWR(cfg_in,SWRfreqs,csc)
%AMSWR Detect sharp wave-ripple events using the discrete Fourier transform.
%
% function SWR = amSWR(cfg_in,SWRfreqs,csc)
%
%   INPUTS
%   SWRfreqs - output from SWRfreak(), noise-corrected frequency spectrum.
%   csc      - output from LoadCSC, comprising the data from a continuously 
%              sampled channel (csc). 
%
%   CONFIGS (with defaults):
%   cfg.weightby = SWRfreqs.parameters.weightby; % 'power' or 'amplitude' 
%              uses the same option that was specified for SWRfreak, but 
%              user can switch 
%              'power' re-weight the spectrum to "unbias" the voltage over 
%               each frequency 
%              'amplitude' use raw spectrum
%   cfg.verbose = 1; % talk to me
%
%   OUTPUTS
%   SWR - tsd with fields:
%         .tvec        - time vector
%         .data        - detection score (similarity) vector
%         .label       - csc used
%         .cfg         - record of cfg history
%         .parameters  - all parameters used to generate the output
%  
% Elyot Grant, Jan 2015 (detection)
% ACarey, Jan 2015 (function design)

%% declare global variables

global cfg parameters
parameters = SWRfreqs.parameters; %because fast hack

%% Parse cfg parameters

cfg_def.verbose = 1;
cfg_def.weightby = SWRfreqs.parameters.weightby;
cfg = ProcessConfig2(cfg_def,cfg_in);

% these settings should be identical to SWRfreak settings 
cfg.hiPassCutoff = SWRfreqs.parameters.hiPassCutoff; % frequency in Hz; delete all freqs below
cfg.fs = SWRfreqs.parameters.fs;

%% Tell the user you're doing something


if cfg.verbose
    tic
    cprintf(-[0 0 1],'amSWR: Identifying potential sharp wave-ripple events...');
    disp(' ');
end

 %% Set up 
 
 %freqs = [0,0,0,0,0,0,0,0,0,0,0,0,0.000170521028425372,0.00170778248757081,0.00331395105024794,0.00479391795042840,0.00607136317217047,0.00678590934283586,0.00697574495989578,0.00650036253423614,0.00569838612070232,0.00466283089230949,0.00366655000685789,0.00272766671935131,0.00195939311675901,0.00137725227252137,0.00102976121026613,0.000775147516298750,0.000613151554386769,0.000445041806482364,0.000371607082692109,0.000310870018623868,0.000321967274810346,0.000312497164011530,0.000308075500895654,0.000273171891686851,0.000204040282703605,0.000109334819915455,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
 
    function swrscore = SWRhelper(freqs,win,csc)

        if strcmp(cfg.weightby,'amplitude') && strcmp(parameters.weightby,'power') % then SWRfreqs were power weighted, but user wants ampl weighted
            freqs = freqs./(1:length(freqs)); 
            
        elseif strcmp(cfg.weightby,'power') && strcmp(parameters.weightby,'amplitude') % then SWRfreqs were ampl weighted, but user wants power weighted
            freqs = freqs.*(1:length(freqs)); 
        end
        
        sampwin = win*cfg.fs; % the window size in nSamples = win * sampling frequency
 
        %hiPassCutoff = 100; %We want to delete all frequencies below 100 Hz
        hiPassCutoff = cfg.hiPassCutoff;
        fourierCoeffCutoff = round(hiPassCutoff*win);
 
        swrscore = zeros(length(csc.data),1);
 
        % FFT
        %for iSWR = sampwin:500000-sampwin % for testing
        for iSWR = sampwin:length(csc.data)-sampwin
        %for iSWR = 3868190:3868191
    
            % Grab a windowed candidate and smooth the ends to reduce noise.
            next = windowedFFT(cfg,csc.data,sampwin,iSWR);
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

if ~isempty(SWRfreqs.freqs2)
    swrscore1 = SWRhelper(SWRfreqs.freqs1,SWRfreqs.parameters.win1,csc);
    % rescale because values are insanely low
    swrscore1 = rescmean(swrscore1,1);
    
    if cfg.verbose
        disp('*yawn* ...') % checkpoint ~ halfway done
    end
    
    swrscore2 = SWRhelper(SWRfreqs.freqs2,SWRfreqs.parameters.win2,csc);
    % rescale
    swrscore2 = rescmean(swrscore2,1);
    
    % combine the scores
    geometricmean = sqrt(swrscore1.*swrscore2);
    
     SWR = tsd(csc.tvec,geometricmean);
else
    swrscore1 = SWRhelper(SWRfreqs.freqs1,SWRfreqs.parameters.win1,csc);
    % rescale because values are insanely low
    swrscore1 = rescmean(swrscore1,1);
    SWR = tsd(csc.tvec,swrscore1);
    
    swrscore2 = [];
end
 
 %% Return output
 mfun = mfilename;
 
 SWR.parameters = struct('weightby',cfg.weightby,'win1',SWRfreqs.parameters.win1,'win2',SWRfreqs.parameters.win2,'hiPassCutoff',cfg.hiPassCutoff,'fs',cfg.fs,'csc',csc.label,'SWRfreak',SWRfreqs.parameters);
 SWR.label = csc.label;
 
 swr1 = tsd(csc.tvec,swrscore1);
 swr1.parameters = struct('weightby',cfg.weightby,'win1',SWRfreqs.parameters.win1,'hiPassCutoff',cfg.hiPassCutoff,'fs',cfg.fs,'csc',csc.label,'SWRfreak',SWRfreqs.parameters);
 swr1.label = csc.label;
 
 swr2 = tsd(csc.tvec,swrscore2);
 swr2.parameters = struct('weightby',cfg.weightby,'win2',SWRfreqs.parameters.win2,'hiPassCutoff',cfg.hiPassCutoff,'fs',cfg.fs,'csc',csc.label,'SWRfreak',SWRfreqs.parameters);
 swr2.label = csc.label;
 
 % keep a record
 if isfield(cfg_in,'history')
     SWR.cfg.history.mfun = cat(1,SWR.cfg.history.mfun,mfilename);
     SWR.cfg.history.cfg = cat(1,SWR.cfg.history.cfg,{cfg});
 else
     SWR.cfg.history.mfun = mfun;
     SWR.cfg.history.cfg = {cfg};
 end
 
 if cfg.verbose
 toc
 end
end