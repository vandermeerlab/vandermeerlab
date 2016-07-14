function csd_out = AMPX_get_CSD_5pt(cfg_in, cycles)
%% AMPX_get_CSD_5pt: computes the csd using the second spatial derivative across sites.
%
%  Uses the csd5pt code created by E.W.Schomburg, March 2014.
%
%
%          Inputs:
%           - cfg_in
%           - cycles [struct] output from the AMPX_Naris_pipeline.  cycles
%           will inlcude the middle/highest amplitude cycle from each gamma
%           event as well as the cycle before and after the it.
%
%          Outputs:
%           - csd_out
%

%% set up parameters
cfg_def.chan_to_process = diag(reshape(1:64, 8,8));
cfg_def.cycle = {'center'};
cfg_def.elec_sep_mm = sqrt(0.2^2+0.2^2); % spacing for A8x8 probe diagonal


cfg = ProcessConfig(cfg_def, cfg_in)


%% loop across all cycles
bands = {'lg', 'hg'};

for iBand = 1:length(bands);
    for iCycle = 1:length(cfg.cycle)
        for iEvt = 1:length(cycles.(bands{iBand}).cycles.peaks.(cfg.cycle{iCycle}).data)
            if isempty(cycles.(bands{iBand}).cycles.peaks.(cfg.cycle{iCycle}).data{iEvt})
                csd_out.(bands{iBand}).(cfg.cycle{iCycle}).csd(:,:,iEvt) = NaN * zeros(4,72); % fill with NaN if the cycle was not deemed usable from the AMPX_get_cycle
            else
                xcpots = cycles.(bands{iBand}).cycles.peaks.(cfg.cycle{iCycle}).data{iEvt}(cfg.chan_to_process,:);
                [csd_out.(bands{iBand}).(cfg.cycle{iCycle}).csd(:,:,iEvt), csd_out.(bands{iBand}).(cfg.cycle{iCycle}).csdelecIds] =csd5pt(xcpots, cfg.elec_sep_mm);
            end
        end
    end
end

%% same for the inverse of the 5th diagonal channel. for relative color axis
for iBand = 1:length(bands);
    for iCycle = 1:length(cfg.cycle)
        for iEvt = 1:length(cycles.(bands{iBand}).cycles.peaks.(cfg.cycle{iCycle}).data)
            if isempty(cycles.(bands{iBand}).cycles.peaks.(cfg.cycle{iCycle}).data{iEvt})
                csd_out.(bands{iBand}).(cfg.cycle{iCycle}).csd(:,:,iEvt) = NaN * zeros(4,72); % fill with NaN if the cycle was not deemed usable from the AMPX_get_cycle
            else
                xcpots = cycles.(bands{iBand}).cycles.peaks.(cfg.cycle{iCycle}).data{iEvt}(cfg.chan_to_process,:);
                xcpots(5,:) = -xcpots(5,:);
                [csd_out.(bands{iBand}).(cfg.cycle{iCycle}).csd_inv(:,:,iEvt), csd_out.(bands{iBand}).(cfg.cycle{iCycle}).csdelecIds_inv] =csd5pt(xcpots, cfg.elec_sep_mm);
            end
        end
    end
end

max_val.lg = max(max(nanmean(csd_out.lg.center.csd_inv, 3)));
max_val.hg = max(max(nanmean(csd_out.hg.center.csd_inv, 3)));

min_val.lg = min(min(nanmean(csd_out.lg.center.csd_inv, 3)));
min_val.hg = min(min(nanmean(csd_out.hg.center.csd_inv, 3)));


%% have a look
figure(100)
subplot(5,1,1:2)
tvec = linspace(0, (cycles.lg.cycles.tvec{1}(end) - cycles.lg.cycles.tvec{1}(1)), length(nanmean(csd_out.lg.center.csd)));
imagesc(tvec,csd_out.lg.center.csdelecIds, nanmean(csd_out.lg.center.csd, 3)); 
set(gca, 'clim', [min_val.lg max_val.lg]);
colorbar('northoutside')

subplot(5,1,3)
plot(nanmean(cell2mat(cycles.lg.cycles.peaks.center.data')))
xlim([1 length(nanmean(cell2mat(cycles.lg.cycles.peaks.center.data')))])

subplot(5,1,4:5)
imagesc(tvec,csd_out.hg.center.csdelecIds, nanmean(csd_out.hg.center.csd, 3));
% set(gca, 'clim', [min_val.hg max_val.hg]);
colorbar('southoutside')




