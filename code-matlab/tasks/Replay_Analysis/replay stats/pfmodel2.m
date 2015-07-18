function [pfm_out] = pfmodel2(cfg_in,pf,pos_in)
%% PFMODEL Place Field Model
%
% Calculates the firing rates for a given set of place cells and time indicies
% based on Chadwick et al., Independent Theta Phase Coding Accounts for CA1 Population
% Sequences and Enables Flexible Remapping. bioRxiv May 12 2014
%
%   INPUTS:
%       cfg_in: cfg file (input parameters)
%       xc: scalar N - creates a single place cell with mean centered at N.
%           1xN array - creates N number of place cells with mean centered at xc(i) values.
%       vel_in: velocity profile tsd
%       
%   CFG OPTIONS:
%       cfg.variance - default 9 (cm).
%           Width of the place field (variance of gaussian rate function r(x)).
%
%       cfg.radius - default 18.75 (cm).
%           Radius of each place field.
%
%       cfg.lockingParameter - default 0.7.
%           Phase locking parameter k. As k->0, the rate tends to r(x). As k->infinity, the
%           rate tends towards r(phi). Note: values above 2 essentially lock the phase to lfp theta          
%
%       cfg.freqTheta - default 8 (Hz).
%           LFP theta value.
%
%       cfg.initTheta - default 0.
%           Initial LFP theta phase offset.
%
%       cfg.deltaPhi - default 2*pi.
%           Total phase precessed over the place field.
%
%       cfg.avgNumSpikes - default 15.
%           Average number of expected spikes for a given place cell per pass of the place
%           field.
%
%       cfg.timebin - default 0.001 (s).
%           The size of timebins. Equivalent to dt.
%
%       cfg.initPos - default 0 cm
%           The initial position of the rat.
%
%       cfg.nTrials - default 100
%           Number of trials to run poisson spike generation
%
%       cfg.tracktype - default "straight"
%           Chooses type of experiment.
%               straight - rat runs in one direction for entire duration
%               linear - rat runs back and forth on a linear track
%
%       cfg.convertfile - default [] returns pfmodel object 
%           Chooses output format.
%               position - position tsd object
%               velocity - velocity tsd object
%               spikes - spike ts object
%               tuning curves - tuning curve tc object
%
%   OUTPUTS:
%       Returns a struct containing the experimental values, a pf (place field)
%       struct and ExpKeys.
%       
%       pf contains parameters for each place cell and the firing rate profile
%       along time vector timestamps. It also contains nTrials generated spiketrains and
%       metadata for these spiketrains. 
%
%       If convertfile is given as input, the output is changed to a specified mvdmlab 
%       formatted object (i.e., ts, tsd)
%
%   EXAMPLE SYNTAX:
%       cfg = [];
%       pfm_out = pfmodel(cfg, [-5 0 5], vel_tsd);
%           Get firing rates for three place cells with means at -5,0,5 cms for a given
%           velocity profile vel_tsd.
%
%
% youkitan 2014-11-02

%% Set defaults and parse optional parameters
cfg_default.variance = 9^2; 
cfg_default.radius = 37.5/2;
cfg_default.lockingParameter = .7; 
cfg_default.freqTheta = 8;
cfg_default.initTheta = 0; 
cfg_default.initPhi = 0; % pf entry at peak of theta
cfg_default.deltaPhi = pi; % pf precess to trough of theta
cfg_default.maxFiringRate = 10; % makes sense if R = sigma
cfg_default.initPos= 0;
cfg_default.nTrials = 100;
cfg_default.speedThreshold = 3; % no spikes when rat is moving less than 3 cm/s
cfg_default.rate_dt = 0.01;
cfg_default.convertfile = [];
cfg_default.verbose = 1;

cfg = ProcessConfig2(cfg_default,cfg_in);
mfun = mfilename;


% Check to see if position tsd has equal sized vectors. Since tvec
% can be different lengths depending on the timebin size, interpolate values for our
% position vector.
assert(all(size(pos_in.data) == size(pos_in.tvec)),...
    'Bad tsd construction: length of vectors do not match!')

tvec = pos_in.tvec;
pvec = pos_in.data;

cfg_spd.window = 0.01;
cfg_spd.postsmoothing = 0.005;
spdvec = abs(dxdt2(cfg_spd,tvec, pvec));

%% Create place cells
for i = 1:length(pf)
%    pf(i).xc = xc(i); 
   pf(i).radius = cfg.radius;
   pf(i).sigma_sq = cfg.variance;
   pf(i).k = cfg.lockingParameter;
   pf(i).initPhi = cfg.initPhi;
   pf(i).deltaPhi = cfg.deltaPhi;
%    pf(i).maxFiringRate = cfg.maxFiringRate;
end

%% Create LFP Theta
initTheta = cfg.initTheta;
freqTheta = cfg.freqTheta;

    function [lfp] = lfpwave(t)
        lfp = cos(initTheta + (2*pi*freqTheta.*t));
    end

%% Create firing rate function
%  Combination of Eq. A1.1 and A2.2 from Chadwick et al. 2014.
%   Normalized so that the maximum value should be A

ratefunc = @(x, theta, pf, v, vthresh) ...
   pf.maxFiringRate ./ exp(pf.k) .* ... % scaled so that the peak value is A
       exp(-(x-pf.xc).^2./(2*pf.sigma_sq))  .* ... % gaussian place field
       exp (pf.k .* cos (pf.initPhi - pf.deltaPhi ./ (2*pf.radius) .* (x - pf.xc + pf.radius) - theta)) .* ... % phase precession
       double(v >= vthresh);

rate_tvec = tvec(1):cfg.rate_dt:tvec(end);
interp_pvec = interp1(tvec,pvec,rate_tvec);
interp_spdvec = interp1(tvec,spdvec,rate_tvec);

if cfg.verbose
    sprintf('\n Creating firing rate profiles...')
end

for nPF=1:length(pf)
    if cfg.verbose; sprintf('.'); end
    
    pf(nPF).rates = ratefunc(interp_pvec, lfpwave(rate_tvec), pf(nPF), ...
      interp_spdvec, cfg.speedThreshold);
end %iterate place cells


%% Generate poisson spikes for N trials
trials = cfg.nTrials;
if cfg.verbose
    sprintf('\n Generate poisson spikes for place cells...')
end
for nPF = 1:length(pf)
    if cfg.verbose; sprintf('.'); end

    pf_func = @(t) ratefunc(interp1(tvec-tvec(1), pvec, t), lfpwave(t), pf(nPF), ...
                      interp1(tvec-tvec(1), spdvec, t), cfg.speedThreshold);

    spiketimes = genInhomogeneousPoisson(cfg.maxFiringRate, tvec(end)-tvec(1), pf_func, trials);
    for N = 1:trials
        spiketimes{N} = spiketimes{N} + tvec(1);
    if isempty(spiketimes{N})
        spiketimes{N} = zeros(0,1);
    end
end %iterate trials

pf(nPF).spiketimes = spiketimes;
end %iterate place cells

%% Housekeeping

% Store place cells and experimental data
data.pf = pf;
data.pvec = pvec;
data.tvec = tvec;
data.vvec = spdvec;
data.rate_tvec = rate_tvec;
data.lfp = lfpwave(tvec);

% Store experimental keys in a single struct
data.ExpKeys.nTrials = trials;
data.ExpKeys.TimeOnTrack = tvec(1);
data.ExpKeys.TimeOffTrack = tvec(end);
data.ExpKeys.TimeBinSize = tvec(2)-tvec(1);

%% Choose output format (can convert to mvdm format)

if ~isempty(cfg.convertfile)
    cases = {'position','velocity','spikes','tuning curves'};
    assert(any(strcmp(cfg.convertfile,cases)),...
        '"%s" is an invalid file type for conversion!',cfg.convertfile)
    switch cfg.convertfile
        case 'position'
            % Create position tsd object
            pos_tsd = tsd;
            pos_tsd.tvec = data.tvec;
            pos_tsd.data = data.pvec;
            pos_tsd.label = {'position (cm)'};
            pos_tsd.cfg.ExpKeys = data.ExpKeys;
            pos_tsd.cfg.history.mfun = cat(1,pos_tsd.cfg.history.mfun,mfun);
            pos_tsd.cfg.history.cfg = cat(1,pos_tsd.cfg.history.cfg,{cfg});
            output = pos_tsd;

        case 'velocity'
            % Create velocity tsd object
            vel_tsd = tsd;
            vel_tsd.tvec = data.tvec;
            vel_tsd.data = data.pvec;
            vel_tsd.label = {'speed (cm/s)'};
            vel_tsd.cfg.ExpKeys = data.ExpKeys;
            vel_tsd.cfg.history.mfun = cat(1,vel_tsd.cfg.history.mfun,mfun);
            vel_tsd.cfg.history.cfg = cat(1,vel_tsd.cfg.history.cfg,{cfg});
            output = vel_tsd;

        case 'spikes'
            % Create spikes ts object
            S = ts;
            for iCell = 1:length(data.pf)
                for iT = 1:cfg.nTrials
                    S.t{iT,iCell} = data.pf(iCell).spiketimes{iT}; %can be multiple trials
                    S.label{iCell} = data.pf(iCell).xc;
                end
            end
            S.cfg.ExpKeys = data.ExpKeys;
            S.cfg.history.mfun = cat(1,S.cfg.history.mfun,mfun);
            S.cfg.history.cfg = cat(1,S.cfg.history.cfg,{cfg});
            output = S;

        case 'tuning curves'
            % Create tuning curve tc object 
            % (Currently not a real tuning curve, just the firing rate profile -- for testing)
            tc_out = tc;
            for iCell = 1:length(data.pf)
                tc_out.firing_rates{iCell} = data.pf(iCell).rates;
                tc_out.stimulus{iCell} = data.pvec;
                tc_out.parameters{iCell} = data.pf(iCell).xc;
                tc_out.label{iCell} = iCell;
            end
            tc_out.cfg.ExpKeys = data.ExpKeys;
            tc_out.cfg.history.mfun = cat(1,tc_out.cfg.history.mfun,mfun);
            tc_out.cfg.history.cfg = cat(1,tc_out.cfg.history.cfg,{cfg});
            output = tc_out;
            
    end %choose output
    pfm_out = output;
    
else
    % Output pfmodel struct
    pfm_out = data;
end





end