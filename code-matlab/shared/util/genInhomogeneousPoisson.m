function [spiketrain] = genInhomogeneousPoisson(max_rate, total_time, ratefunc, nTrials);
% function [spiketrain] = genInhomogeneousPoisson(max_rate, total_time, ratefunc);
%  Generate an inhomogenous Poisson process spike train using deletion
%   max_rate = maximum instantaneous firing rate
%   total_time = length of time to simulate for (starts at time 0)
%   ratefunc = function which should return rate at time t (should not exceed max_rate)

% Generate inhomogeneous Poisson process

if nargin < 3
  error('Must define rate function.');
elseif nargin < 4
  nTrials = 1;
end

Homogeneous_Spikes = poissrnd( (max_rate * total_time) * ones(nTrials,1) );
% beta = 1.2;
% Homogeneous_Spikes = round(gamrnd((max_rate * total_time)/beta, beta) * ones(nTrials,1) );

for k = 1:nTrials
    Homogeneous_Spike_Times = rand(Homogeneous_Spikes(k),1) * total_time;
    spikes = sort(Homogeneous_Spike_Times);
    thinning_prob = rand(length(spikes),1);
    relative_intensity = (1/max_rate) * ratefunc(spikes);
    deleted_spikes = thinning_prob > relative_intensity;
    spikes(deleted_spikes) = [];
    spiketrain{k} = spikes;
%     spiketrain{k} = length(spikes);
end