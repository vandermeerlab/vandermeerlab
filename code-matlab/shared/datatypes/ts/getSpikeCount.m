function spk_c = getSpikeCount(cfg_in,S)
% function spk_c = getSpikeCount(cfg_in,S)
%
% returns spike counts (number of spikes) in a nCells x 1 vector of S input
% 
% if optional cfg_in.iv argument is specified, return spike counts for each
% interval (output is nCells x nIntervals)
%
% MvdM 2016-01-11 initial version

cfg_def = [];
cfg_def.iv = []; % if specified, return spike counts for each interval (output is nCells x nIntervals)
cfg_def.convertToRate = 0; % if 1, divide spike count by length of interval
cfg_def.verbose = 1;

cfg = ProcessConfig(cfg_def,cfg_in);

if ~CheckTS(S)
   error('Input is not a correctly formed ts.'); 
end

if isempty(cfg.iv) % create dummy iv to avoid repetition later
   cfg.iv = iv();
   cfg.iv.tstart = firstSpike(S)-eps;
   cfg.iv.tend = lastSpike(S)+eps;
end

nCells = length(S.t);

for iI = length(cfg.iv.tstart):-1:1 % for each interval, get spike counts
    
   this_S = restrict(S,cfg.iv.tstart(iI),cfg.iv.tend(iI));
   
   for iC = nCells:-1:1
  
       spk_c(iC,iI) = length(this_S.t{iC});
       
       if cfg.convertToRate, spk_c(iC,iI) = spk_c(iC,iI) ./ (cfg.iv.tend(iI) - cfg.iv.tstart(iI)); end
       
   end  % of cells
    
end % of intervals