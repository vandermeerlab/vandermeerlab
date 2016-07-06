function [h1, h2] = AMPX_Naris_fig_3_phase(cfg_in, pre, post, varargin)
%%AMPX_Naris_fig_3_phase: creates two 4 x 2 figues of the phase differences
% for both across entire gamma events and across the central three cycles.
%
%          Inputs:
%           - cfg_in: [struct] contains parameters
%           - all_data: [struct] all the data from each session (output
%           from AMPX_Naris_pipeline).
%           -
%          Outputs:
%           - h1: handle for the across events figure
%           - h2: handle for the across triplet figure
%
% EC - 2016-05-23

%% set up parameters
extract_varargin

cfg_def = [];

cfg = ProcessConfig2(cfg_in, cfg_def);


%% get the average phase across sessions for each rat
% initialize
sess_list = fieldnames(pre);
bands = {'lg', 'hg'};
for iBand = 1:length(bands)
    R061.(bands{iBand}).phase = []; R061.(bands{iBand}).cycles = [];
    R054.(bands{iBand}).phase = []; R054.(bands{iBand}).cycles = [];
    R045.(bands{iBand}).phase = []; R045.(bands{iBand}).cycles = [];
    R049.(bands{iBand}).phase = []; R049.(bands{iBand}).cycles = [];
end

for iSess = 1:length(sess_list)
    for iBand = 1:length(bands)
        if strcmp(sess_list{iSess}(1:4), 'R061')
            R061.(bands{iBand}).phase = cat(3,R061.(bands{iBand}).phase, pre.(sess_list{iSess}).(bands{iBand}).phase.phase_avg_mean);
            R061.(bands{iBand}).cycles = cat(3,R061.(bands{iBand}).cycles, pre.(sess_list{iSess}).(bands{iBand}).cycles_phase.phase.avg.mean);
        end
        if strcmp(sess_list{iSess}(1:4), 'R054')
            R054.(bands{iBand}).phase = cat(3,R054.(bands{iBand}).phase, pre.(sess_list{iSess}).(bands{iBand}).phase.phase_avg_mean);
            R054.(bands{iBand}).cycles = cat(3, R054.(bands{iBand}).cycles, pre.(sess_list{iSess}).(bands{iBand}).cycles_phase.phase.avg.mean);
        end
        if strcmp(sess_list{iSess}(1:4), 'R045')
            R045.(bands{iBand}).phase = cat(3,R045.(bands{iBand}).phase, pre.(sess_list{iSess}).(bands{iBand}).phase.phase_avg_mean,post.(sess_list{iSess}).(bands{iBand}).phase.phase_avg_mean );
            R045.(bands{iBand}).cycles = cat(3,R045.(bands{iBand}).cycles, pre.(sess_list{iSess}).(bands{iBand}).cycles_phase.phase.avg.mean, post.(sess_list{iSess}).(bands{iBand}).cycles_phase.phase.avg.mean);
        end
        if strcmp(sess_list{iSess}(1:4), 'R049')
            R049.(bands{iBand}).phase = cat(3,R049.(bands{iBand}).phase, pre.(sess_list{iSess}).(bands{iBand}).phase.phase_avg_mean, post.(sess_list{iSess}).(bands{iBand}).phase.phase_avg_mean);
            R049.(bands{iBand}).cycles = cat(3,R049.(bands{iBand}).cycles, pre.(sess_list{iSess}).(bands{iBand}).cycles_phase.phase.avg.mean, post.(sess_list{iSess}).(bands{iBand}).cycles_phase.phase.avg.mean);
        end
    end
end


%% average acros rats

R54_phase.lg = rad2deg(mean(R054.lg.phase, 3)); R54_phase.hg = rad2deg(mean(R054.hg.phase, 3));
R54_cycles.lg = rad2deg(mean(R054.lg.cycles, 3)); R54_cycles.hg = rad2deg(mean(R054.hg.cycles, 3));

R61_phase.lg = rad2deg(mean(R061.lg.phase, 3)); R61_phase.hg = rad2deg(mean(R061.hg.phase, 3));
R61_cycles.lg = rad2deg(mean(R061.lg.cycles, 3)); R61_cycles.hg = rad2deg(mean(R061.hg.cycles, 3));

R45_phase.lg = rad2deg(mean(R045.lg.phase, 3)); R45_phase.hg = rad2deg(mean(R045.hg.phase, 3));
R45_cycles.lg = rad2deg(mean(R045.lg.cycles, 3)); R45_cycles.hg = rad2deg(mean(R045.hg.cycles, 3));

R49_phase.lg = rad2deg(mean(R049.lg.phase, 3)); R49_phase.hg = rad2deg(mean(R049.hg.phase, 3));
R49_cycles.lg = rad2deg(mean(R049.lg.cycles, 3)); R49_cycles.hg = rad2deg(mean(R049.hg.cycles, 3));

%% get the min and max for the whole set, used for Clims in the imagesc
lg_min = nanmin(nanmin([R54_phase.lg, R61_phase.lg, R45_phase.lg, R49_phase.lg]));
hg_min = nanmin(nanmin([R54_phase.hg, R61_phase.hg, R45_phase.hg, R49_phase.hg]));

lg_min_cycle = nanmin(nanmin([R54_cycles.lg, R61_cycles.lg, R45_cycles.lg, R49_cycles.lg]));
hg_min_cycle = nanmin(nanmin([R54_cycles.hg, R61_cycles.hg, R45_cycles.hg, R49_cycles.hg]));

lg_max = nanmax(nanmax([R54_phase.lg, R61_phase.lg, R45_phase.lg, R49_phase.lg]));
hg_max = nanmax(nanmax([R54_phase.hg, R61_phase.hg, R45_phase.hg, R49_phase.hg]));

lg_max_cycle = nanmax(nanmax([R54_cycles.lg, R61_cycles.lg, R45_cycles.lg, R49_cycles.lg]));
hg_max_cycle = nanmax(nanmax([R54_cycles.hg, R61_cycles.hg, R45_cycles.hg, R49_cycles.hg]));
%% make the figure
% across events
h1 = figure(11);
set(h1, 'PaperPositionMode', 'auto', 'color', 'w')
set(h1,  'position', [263, -19, 1307, 585])
% low gamma on top, high on bottom, rat order 61, 54, 45, 49
subplot(2,4,1)
nan_imagesc_ec(R61_phase.lg);
subplot(2,4,2)
nan_imagesc_ec(R54_phase.lg);
subplot(2,4,3)
nan_imagesc_ec(R49_phase.lg);
subplot(2,4,4)
nan_imagesc_ec(R45_phase.lg);

subplot(2,4,5)
nan_imagesc_ec(R61_phase.hg);
subplot(2,4,6)
nan_imagesc_ec(R54_phase.hg);
subplot(2,4,7)
nan_imagesc_ec(R49_phase.hg);
subplot(2,4,8)
nan_imagesc_ec(R45_phase.hg);

for ii = 1:8
    subplot(2,4,ii)
    set(gca, 'xtick', [], 'ytick', []);
    if ii < 5
        set(gca, 'clim', [-10 10])
    else
        set(gca, 'clim', [-10 10])
    end
end
% hp4 = get(subplot(2,4,8),'Position');
% colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])
%% across cycles
h2 = figure(22);
set(h2, 'PaperPositionMode', 'auto', 'color', 'w')
set(h2,  'position', [263, -19, 1307, 585])
% low gamma on top, high on bottom, rat order 61, 54, 45, 49
subplot(2,4,1)
nan_imagesc_ec(R61_cycles.lg)
subplot(2,4,2)
nan_imagesc_ec(R54_cycles.lg)
subplot(2,4,3)
nan_imagesc_ec(R49_cycles.lg)
subplot(2,4,4)
nan_imagesc_ec(R45_cycles.lg)

subplot(2,4,5)
nan_imagesc_ec(R61_cycles.hg)
subplot(2,4,6)
nan_imagesc_ec(R54_cycles.hg)
subplot(2,4,7)
nan_imagesc_ec(R49_cycles.hg)
subplot(2,4,8)
nan_imagesc_ec(R45_cycles.hg)

for ii = 1:8
    subplot(2,4,ii)
    set(gca, 'xtick', [], 'ytick', []);
    %     colorbar
    if ii < 5
        set(gca, 'clim', [-10 10])
    else
        set(gca, 'clim', [-10 10])
    end
end


%% export the figures
if exist('save_fig')
    saveas(h1, 'D:\DATA\Paper_figs\Fig3_A', 'fig')
    saveas(h1, 'D:\DATA\Paper_figs\Fig3_A', 'epsc')
    
    saveas(h2, 'D:\DATA\Paper_figs\Fig3_B', 'fig')
    saveas(h2, 'D:\DATA\Paper_figs\Fig3_B', 'epsc')
end
