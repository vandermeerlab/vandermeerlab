function S = RemoveInterneuronsHC(cfg_in,S,lfp)
% function S = RemoveInterneuronsHC(cfg_in,S,lfp)
%
% removes neurons with mean firing rate above cfg.max_fr
%
% MvdM 2016

cfg_def = [];
cfg_def.showFRhist = 0;
cfg_def.max_fr = 5;

cfg = ProcessConfig(cfg_def,cfg_in);

% use LFP to calculate true experiment time -- need this because of
% possible gaps in recording, so can't just take first and last spike
lfp_dt = median(diff(lfp.tvec));
total_exp_time = length(lfp.tvec).*lfp_dt;

nCells = length(S.t);
for iC = nCells:-1:1
    
   this_t = S.t{iC};
   this_fr(iC) = length(this_t)./total_exp_time;
    
end

if cfg.showFRhist
   hist(this_fr,100); 
end

keep_idx = this_fr <= cfg.max_fr;

S = SelectTS([],S,keep_idx);

fprintf('RemoveInterneuronsHC: %d/%d neurons kept.\n',sum(keep_idx),nCells);
