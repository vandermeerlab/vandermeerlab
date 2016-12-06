function [spiketrain] = genInhomogeneousPoisson(max_rate, total_time, ratefunc, nTrials)
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

Homogeneous_Spikes = poissrnd( (max_rate * total_time) * ones(nTrials,1) ); %number of spikes
% beta = 1.2;
% Homogeneous_Spikes = round(gamrnd((max_rate * total_time)/beta, beta) * ones(nTrials,1) );

for k = 1:nTrials
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
%     spiketrain{k} = length(spikes);
end

end

function spikes_out = renewal(spikes_in)
    del_idx = find(diff(spikes_in) < 0.001); %1ms refractory
    if ~isempty(del_idx)
        del_idx = del_idx + 1; %this is actual index of successive spike time
        spikes_in(del_idx) = [];
    end
    spikes_out = spikes_in;
end
