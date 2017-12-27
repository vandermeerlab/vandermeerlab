function tsd_out = triclops(cfg_in,tsd_in)
%TRICLOPS Perform local adaptive thresholding on TSD
%  TSD = triclops(cfg,TSD) Compares activity level across three viewing
%    windows and rescales the data point in the center according to the
%    algorithm chosen by the user.
%
%    INPUTS
%       cfg: config struct controlling function behaviour
%       TSD: timestamped data such as MUA
%
%    OUTPUTS:
%       TSD: timestamped data rescaled using adaptive thresholding
%
%    CONFIG OPTIONS
%       cfg.Cwin = 120; Width of the central window, in *nSamples*
%       cfg.Nwin = 720; Width of the neighbouring windows in *nSamples*
%       cfg.method = 'percerr'
%       cfg.stepSize = 4; How often to compute the score, once every n
%                samples for indirect speed control. Intermediate data
%                points are interpolated.
%       cfg.InterpMethod = 'pchip'; The method for interpolating points
%                skipped as a result of cfg.stepSize.
%       cfg.verbose = 1; If 1 function prints text to command window, if 0
%                        is silent.
%
% aacarey Aug 2015, edit Jan 2016
%
% see also: amMUA

%% Parse cfg parameters
tic
cfg_def.Cwin = 120;
cfg_def.Nwin = 720;
cfg_def.method = 'percerr';
cfg_def.InterpMethod = 'pchip';
cfg_def.stepSize = 4;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

disp([mfun,': performing adaptive rescaling...'])

%% Do the thing
% we want to compare regions to their neighbours and generate a new vector
% that describes how active a localized region is in comparison to its
% near neighbours. This allows us to do adaptive thresholding. If it's done
% well engough, we should be able to detect SWR-associated spiking activity
% when it occurs amidst theta-associated spiking activity.

% one intuitive way of comparing a region to its neighbours would be to use
% a central viewing window (for the region of interest) and two surrounding
% windows (for its immediate neighbours). Depending on the relative
% activity accross the three regions, we assign the focal region a score. 

%Cwin = 120; % center window in nSamples. If fs = 2000, then this is a 60 ms time window
Cwin = cfg.Cwin;
%Nwin = 720; % neighbour sample window
Nwin = cfg.Nwin;
halfspan = (Cwin + 2*Nwin)/2; % consider all three windows, halfspan is the distance in nSamples from the center of the
% center viewing window out to the edge of the neighbouring windows.
% halfspan is necessary so we know where in the vector to start and stop
% scoring.

score = NaN(length(tsd_in.tvec),1); % preallocate space for speed
switch cfg.method
    case 'ratio'
        % ratio: score = c/left + c/right
        for iFocus = halfspan+1:cfg.stepSize:length(tsd_in.tvec)-halfspan
            
            % construct each of triclops' eyes
            Leye = [iFocus-halfspan iFocus-Cwin/2];
            Ceye = [iFocus-Cwin/2 iFocus+Cwin/2];
            Reye = [iFocus+Cwin/2 iFocus+halfspan];
            
            % grab the MUA that triclops sees. the activity is the total number
            % of spikes that occur within each eye (viewing window). but because the windows
            % may be different sizes, weight the scores accordingly by dividing by
            % the window width
            
            Lactivity = sum(tsd_in.data(Leye(1):Leye(2)))/Nwin + 1; % +1 to avoid dividing by zero
            Cactivity = sum(tsd_in.data(Ceye(1):Ceye(2)))/Cwin + 1;
            Ractivity = sum(tsd_in.data(Reye(1):Reye(2)))/Nwin + 1;
        
            % generate a score for the focal point, based on the activity of the
            % center region compared to its neighbours. this is like an "adaptive"
            % score
            
            focalscore = Cactivity/Lactivity + Cactivity/Ractivity;
            if isnan(focalscore); oioi; end
            score(iFocus) = focalscore;
        end
    case 'percerr'       
        % calc like percentage error, but the constant 100 is stupid
        % score = (c-left)/left*100 + (c-right)/right*100
        for iFocus = halfspan+1:cfg.stepSize:length(tsd_in.tvec)-halfspan
            % construct each of triclops' eyes
            Leye = [iFocus-halfspan iFocus-Cwin/2];
            Ceye = [iFocus-Cwin/2 iFocus+Cwin/2];
            Reye = [iFocus+Cwin/2 iFocus+halfspan];
            
            Lactivity = sum(tsd_in.data(Leye(1):Leye(2)))/Nwin;
            Cactivity = sum(tsd_in.data(Ceye(1):Ceye(2)))/Cwin;
            Ractivity = sum(tsd_in.data(Reye(1):Reye(2)))/Nwin;
            
            % generate score for the center window
            focalscore = (Cactivity-Lactivity)/(Lactivity+1) + (Cactivity-Ractivity)/(Ractivity+1);
            score(iFocus) = focalscore;
            % i add 1 to avoid the situation where x is a very small number
            % in 1/x; rather i use 1/(1+x). This prevents border regions
            % from going to infinity (areas where there's a spike in the
            % center window but zero spikes in the neighbour)
        end
        
    case 'am1' % not good
        % calc like percentage error, but the constant 100 is stupid
        % score = (c-left)/left*100 + (c-right)/right*100
        for iFocus = halfspan+1:cfg.stepSize:length(tsd_in.tvec)-halfspan
            % construct each of triclops' eyes
            Leye = [iFocus-halfspan iFocus-Cwin/2];
            Ceye = [iFocus-Cwin/2 iFocus+Cwin/2];
            Reye = [iFocus+Cwin/2 iFocus+halfspan];
            
            Lactivity = sum(tsd_in.data(Leye(1):Leye(2)))/Nwin;
            Cactivity = sum(tsd_in.data(Ceye(1):Ceye(2)))/Cwin;
            Ractivity = sum(tsd_in.data(Reye(1):Reye(2)))/Nwin;
            
            % generate score for the center window
            focalscore = Cactivity/(Lactivity+Ractivity+1);
            score(iFocus) = focalscore;
        end
        
    case 'am2' % not good
        % calc like percentage error, but the constant 100 is stupid
        % score = (c-left)/left*100 + (c-right)/right*100
        for iFocus = halfspan+1:cfg.stepSize:length(tsd_in.tvec)-halfspan
            % construct each of triclops' eyes
            Leye = [iFocus-halfspan iFocus-Cwin/2];
            Ceye = [iFocus-Cwin/2 iFocus+Cwin/2];
            Reye = [iFocus+Cwin/2 iFocus+halfspan];
            
            Lactivity = sum(tsd_in.data(Leye(1):Leye(2)))/Nwin;
            Cactivity = sum(tsd_in.data(Ceye(1):Ceye(2)))/Cwin;
            Ractivity = sum(tsd_in.data(Reye(1):Reye(2)))/Nwin;
            
            % generate score for the center window
            focalscore = (2*Cactivity-Lactivity-Ractivity)/(Lactivity+Ractivity+1);
            score(iFocus) = focalscore;
        end
        
    case 'am3' % not good
        for iFocus = halfspan+1:cfg.stepSize:length(tsd_in.tvec)-halfspan
            % construct each of triclops' eyes
            Leye = [iFocus-halfspan iFocus-Cwin/2];
            Ceye = [iFocus-Cwin/2 iFocus+Cwin/2];
            Reye = [iFocus+Cwin/2 iFocus+halfspan];
            
            Lactivity = sum(tsd_in.data(Leye(1):Leye(2)))/Nwin;
            Cactivity = sum(tsd_in.data(Ceye(1):Ceye(2)))/Cwin;
            Ractivity = sum(tsd_in.data(Reye(1):Reye(2)))/Nwin;
            
            % generate score for the center window
            focalscore = 2*Cactivity/(Lactivity+Ractivity+1);
            score(iFocus) = focalscore;
            
        end
    case '!@#$' % censored
        for iFocus = halfspan+1:cfg.stepSize:length(tsd_in.tvec)-halfspan
            % construct each of triclops' eyes
            Leye = [iFocus-halfspan iFocus-Cwin/2];
            Ceye = [iFocus-Cwin/2 iFocus+Cwin/2];
            Reye = [iFocus+Cwin/2 iFocus+halfspan];
            
            Lactivity = sum(tsd_in.data(Leye(1):Leye(2)))/Nwin;
            Cactivity = sum(tsd_in.data(Ceye(1):Ceye(2)))/Cwin;
            Ractivity = sum(tsd_in.data(Reye(1):Reye(2)))/Nwin;
            
            % generate score for the center window
            focalscore = Cactivity;
            score(iFocus) = focalscore;
        end
    case 'zscore'
        % for zscore...no confidence in whether this is correct, but it
        % takes forever to run
        for iFocus = halfspan+1:cfg.stepSize:length(tsd_in.tvec)-halfspan
            Leye = [iFocus-halfspan iFocus-Cwin/2];
            Ceye = [iFocus-Cwin/2 iFocus+Cwin/2];
            Reye = [iFocus+Cwin/2 iFocus+halfspan];
            Lvals = tsd_in.data(Leye(1):Leye(2));
            Cvals = tsd_in.data(Ceye(1):Ceye(2));
            Rvals = tsd_in.data(Reye(1):Reye(2));
            Lactivity = mean(Lvals)/Nwin;
            Cactivity = mean(Cvals)/Cwin;
            Ractivity = mean(Rvals)/Nwin;
            
            % generate a score for center window
            focalscore = (Cactivity-Lactivity)/std(Lvals) + (Cactivity-Ractivity)/std(Rvals); 
            score(iFocus) = focalscore;
        end
        
    otherwise
        error('Unrecognized method. Better check that spelling.')
end

if cfg.stepSize > 1 % then interpolate the intermediate values
    score(1) = 0; score(end) = 0; % because hack or interp can mess up the flanks
    NonNaNs = ~isnan(score); % these are indices of the values returned by the loop above
    
    % now interpolate the missing sample points
    score = interp1(tsd_in.tvec(NonNaNs),score(NonNaNs),tsd_in.tvec,cfg.InterpMethod); % can do linear too
end

%% convert to tsd

tsd_out = tsd(tsd_in.tvec,score');
tsd_out = History(tsd_out,mfilename,cfg);

toc

end

