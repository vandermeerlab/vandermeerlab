function [all_cycles] = AMPX_phase_cycle(cfg_in, all_cycles)
%% AMPX_phase_cycle: uses the extracted gamma peak cycles from the AMPX_kluster_kCSD
% function and uses the cpsd method to extract the phase for the three
% central cycle sin the events.
%
%
%inputs:
%         - cfg: configs.  Should contain window parameters
%         - cycles:  the all_cycles output from AMPX_kluster_kCSD
%
%Outputs
%
%         - all_cycles [struct] adds the phase subfield with the cpsd and
%         cfg.
%
%%%

%% extract the phase differences across events.
cfg_def.debug = 0;
cfg_def.spec_window = 512;
cfg_def.bands = {'low', 'high'};
cfg_def.ref = Naris_BestChan_remap(all_cycles.ExpKeys, 'location', 'vl');
cfg = ProcessConfig2(cfg_def, cfg_in);
cfg.date_processed = datestr(now);



% for iband = 1:length(cfg.bands); % loop over low and then high bands.  Can be just one if speficied in the cfg_in

%     cfg_band = AMPX_which_band(cfg.bands{iband}); % gets the band specific parameters.  Keeps things consistent across fcns


for ievt = length(all_cycles.data):-1:1 % loop over the events
    if isempty(all_cycles.data{ievt})
        all_cycles.phase.cpsd(:,:,ievt) = ones(8,8)*NaN;
        continue
    else
        cfg.fs =  1/mode(diff(all_cycles.tvec{1,ievt}));
        
        
        for iChan = 64:-1:1
            if ismember(all_cycles.ExpKeys.Probe_layout(iChan), all_cycles.ExpKeys.BadChannels)
                coh_specs(iChan) = NaN;
            else
                [Cxy,F] = cpsd(all_cycles.data{ievt}(cfg.ref,:),all_cycles.data{ievt}(iChan,:),hanning(cfg.spec_window/4),cfg.spec_window/8,1024,cfg.fs);
                %     [Cxy,F] = cpsd(cycles(1,:),cycles(iChan,:),hamming(cfg.spec_window/4),cfg.spec_window/8,1024,cfg.fs);
                %     [Cxy,F] = cpsd(cycles(1,:), cycles(iChan,:), hamming(.25*length(cycles(1,:))), .125*length(cycles(1,:)), cfg.fs);
                coh_spec = -angle(Cxy); %higher value means leading. outputs radians
                coh_specs(iChan) = circ_mean(coh_spec(nearest(F,cfg.f(1)):nearest(F,cfg.f(2))));
            end
        end
    end
    all_cycles.phase.cpsd(:,:,ievt) = reshape(rad2deg(coh_specs),8,8);
    temp_for_avg = rad2deg(coh_specs);
    temp_for_avg(isnan(temp_for_avg)) =0;
    all_cycles_for_avg(:,:, ievt) = reshape(temp_for_avg,8,8);
    
end

all_cycles.phase.avg.mean = circ_mean(degtorad(all_cycles_for_avg),[],3);
for x = 1:8
    for y = 1:8
        all_cycles.phase.avg.median(x,y) = circ_median(degtorad(all_cycles_for_avg(x,y,:)));
    end
end

all_cycles.phase.avg.mean(all_cycles.phase.avg.mean ==0) = NaN;
all_cycles.phase.avg.median(all_cycles.phase.avg.median ==0) = NaN;
all_cycles.phase.avg.mean(cfg.ref) = 0;
all_cycles.phase.avg.median(cfg.ref) = 0;
all_cycles.phase.cfg = cfg;
%%
close all
figure(10)
imagescnan(all_cycles.phase.avg.mean)
colorbar('location','eastoutside')
caxis([-(max(max(abs(all_cycles.phase.avg.mean)))), (max(max(abs(all_cycles.phase.avg.mean))))])
set(gca, 'xtick',[], 'ytick',[])
ax = axes('position',[0,0,1,1],'visible','off');
tx = text(0.46,0.95,[all_cycles.ExpKeys.fname '_avg']);
set(tx,'fontweight','bold');set(tx,'fontsize',10);
set(gcf,'Units','pixels','Position',[680   558   560   420])

%% save a figure

%     saveas(gcf, ['C:\Users\mvdmlab\Dropbox\Naris_Paper\Phase_figure\' strrep(all_cycles.ExpKeys.fname, '-', '_') '_avg_phase.fig'])
if cfg.debug
    for ievt = 1:length(all_cycles.phase.cpsd)
        if abs(diff([max(max(all_cycles.phase.cpsd(:,:,ievt))), min(min(all_cycles.phase.cpsd(:,:,ievt)))])) > 20
            imagescnan(all_cycles.phase.cpsd(:,:,ievt))
            all_cycles.phase.cpsd(:,:,ievt)
            disp(num2str(ievt))
            pause
        end
    end
end
%% check
% reshape(all_cycles.ExpKeys.Probe_layout, 8,8)
% reshape(rad2deg(coh_specs),8,8)
% mean(all_cycles.(cfg.bands{iband}).phase.cpsd, 3)