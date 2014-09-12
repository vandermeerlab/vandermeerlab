function Q = MakeQfromS(cfg_in,S)
% function Q = MakeQfromS(cfg,S)
%
% Make Q matrix from ts with spike data
%
% INPUTS:
%
% OUTPUTS:
%
% MvdM 2014-08-21

cfg = [];
cfg.dt = 0.05;
cfg.mode = 'spikeCount'; % 'firingRate'
cfg.tvec_edges = [];

ProcessConfig;
mfun = mfilename;

% assemble tvecs
if isempty(cfg.tvec_edges)
    cfg.tvec_edges = firstSpike(S):cfg.dt:lastSpike(S);
else
    cfg.dt = median(diff(cfg.tvec_edges));
end

tvec_centers = cfg.tvec_edges(1:end-1)+cfg.dt/2;

switch cfg.mode
    case 'spikeCount'
        for iC = length(S.t):-1:1
            
            spk_t = S.t{iC};
            Q(iC,:) = trim_histc(histc(spk_t,cfg.tvec_edges));
            
        end
    case 'firingRate'
        error('firingRate mode not yet implemented.');
end

Q = tsd(tvec_centers,Q);

% housekeeping
Q.cfg.history.mfun = cat(1,Q.cfg.history.mfun,mfun);
Q.cfg.history.cfg = cat(1,Q.cfg.history.cfg,{cfg});
