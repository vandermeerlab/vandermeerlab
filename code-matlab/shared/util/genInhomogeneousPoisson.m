function spiketrain = genInhomogeneousPoisson(cfg_in,max_rate, total_time, ratefunc)
%% GENINHOMOGENEOUSPOISSON Generate Inhomogeneous Poisson spikes 
% [spiketrain] = genInhomogeneousPoisson(max_rate, total_time, ratefunc) Generates an
% inhomogenous Poisson process spike train using a deletion method
%
%   INPUT:
%       max_rate - maximum instantaneous firing rate
%       total_time - length of time to simulate for (starts at time 0)
%       ratefunc - function which should return rate at time t (should not exceed max_rate)
%
%   OUTPUT:
%       S_out - output ts object
%
%   CFG OPTIONS:
%       cfg.nTrials = 1; Number of trials.
%       cfg.rp = .001; Refractory period in seconds.
%       cfg.verbose = 0; If 1, print informative text to command window; if 0, be silent.
%
% youkitan 2016-12-08 edits: cfg parameters

%% parse cfg parameters
cfg_def = [];
cfg_def.nTrials = 1;
cfg_def.refractoryPeriod = .001;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

%% main body

Homogeneous_Spikes = poissrnd( (max_rate * total_time) * ones(cfg.nTrials,1) ); %number of spikes
% beta = 1.2;
% Homogeneous_Spikes = round(gamrnd((max_rate * total_time)/beta, beta) * ones(nTrials,1) );

for k = 1:cfg.nTrials
    Homogeneous_Spike_Times = rand(Homogeneous_Spikes(k),1) * total_time; %convert to spiketimes
    spikes = sort(Homogeneous_Spike_Times);
    spikes = renewal(spikes); % enforced absolute refractory period of 1ms
    thinning_prob = rand(length(spikes),1); %pick unif random threshold
    relative_intensity = (1/max_rate) * ratefunc(spikes); %normalized spike prob
    deleted_spikes = thinning_prob > relative_intensity; %equal to keep = probspike > threshold
    spikes(deleted_spikes) = [];
    
    if isempty(spikes)  %this ensures compatibility with spike plotting and fails
        spikes = zeros(0,1);
    end
    spiketrain{k} = spikes;
end

end

function spikes_out = renewal(spikes_in)
    del_idx = find(diff(spikes_in) < cfg.refractoryPeriod); %1ms refractory
    if ~isempty(del_idx)
        del_idx = del_idx + 1; %this is actual index of successive spike time
        spikes_in(del_idx) = [];
    end
    spikes_out = spikes_in;
end
