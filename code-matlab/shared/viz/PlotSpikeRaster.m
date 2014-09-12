function PlotSpikeRaster(cfg_in,S)
% function PlotSpikeRaster(cfg,S)
%
% MvdM 2014-07-20

cfg.SpikeHeight = 0.4;
ProcessConfig;

for iC = 1:length(S.t)
    
    if isempty(S.t{iC})
        continue;
    end
    
    h{iC} = plot([S.t{iC} S.t{iC}],[iC-cfg.SpikeHeight iC+cfg.SpikeHeight],'k');
    hold on;
    
end


if isfield(cfg,'lfp')
    
   cfg.lfp.data = rescale(cfg.lfp.data,-5,0);
   plot(cfg.lfp.tvec,cfg.lfp.data,'Color',[0.5 0.5 0.5]);
    
end

axis tight;
