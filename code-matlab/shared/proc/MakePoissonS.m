function S = MakePoissonS(cfg_in,sdf)
% function S = MakePoissonS(cfg_in,sdf)
%
% generates Poisson spike trains
%
% INPUTS:
% sdf: 1 x nCells array with mean firing rates
%
% CONFIGS:
% cfg_def.dt = 0.001; % time step when generating spikes (in s)
% cfg_def.t = [0 10]; % interval to generate spikes for (in s)
%
% MvdM 2016-01-18 initial version

cfg_def = [];
cfg_def.dt = 0.001; % time step when generating spikes (in s)
cfg_def.t = [0 10]; % interval to generate spikes for (in s)
% need option for setting RNG to specific state
% need option for different algorithms to generate spikes

cfg = ProcessConfig(cfg_def,cfg_in);
mfun = mfilename;

tvec = cfg.t(1):cfg.dt:cfg.t(2);

% if cdf is a tsd: inhomogenous Poisson process
if CheckTSD(sdf)
    error('Inhomogenous Poisson process not yet implemented.');
else % if cdf is a vector: homogenous Poisson process (one mean FR for each cell)

    S = ts;
    
    for iC = 1:length(sdf)
              
        pspike = sdf(iC)*cfg.dt; % probability of generating a spike in bin

        spk_poiss = rand(size(tvec)); % random numbers between 0 and 1
        spk_poiss_idx = find(spk_poiss < pspike); % index of bins with spike
        
        S.t{iC} = tvec(spk_poiss_idx)'; % use idxs to get corresponding spike time
        
    end
    
end

% housekeeping
S.cfg.history.mfun = cat(1,S.cfg.history.mfun,mfun);
S.cfg.history.cfg = cat(1,S.cfg.history.cfg,{cfg});