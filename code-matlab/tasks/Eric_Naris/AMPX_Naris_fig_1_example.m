function [h1, h_sub2] = AMPX_Naris_fig_1_example(all_data_pre, all_data_post, varargin)
%%AMPX_Naris_fig_3_phase: creates two 4 x 2 figues of the phase differences
% for both across entire gamma events and across the central three cycles.
%
%          Inputs:
%           - cfg_in: [struct] contains parameters
%           - all_data_pre: [struct] all the data from each 'pre' session (output
%           from AMPX_Naris_pipeline).
%           - al_data_post: [struct] all the data from each 'post' session (output
%           from AMPX_Naris_pipeline).
%          Outputs:
%           - h1: handle for the example raw trace/heat map  Fig1 C/D
%           - h2: handle for the average across rats  Fig1 E
%           - h3: handle for the average across rats (min/max c-axis) S1
%
% EC - 2016-05-23

%% set up parameters

cfg_def = [];
cfg_def.markers = {'#', '+', 'x', 'o'};
cfg_def.mrk_off = -5;
cfg_def.raw_plot_boost = 300;
cfg_def.line_end  =.5;
cfg_def.width = .5;
cfg_def.ft_size = 20;
cfg_def.linewidth = 2;
cfg_def.example = [];
cfg_def.session_name = [];
cfg = ProcessConfig(cfg_def, cfg_in);

if isempty(cfg.session_name) ~=1
    cfg.session_name = 'R061_2014_09_26'; % most channels & best positioning
end
%% load the data file to get the raw traces
fname = strrep(cfg.session_name, '_', '-');
cd(['D:\DATA\' fname(1:4) '\' fname(1:15) ])
fname = strrep(cfg.session_name, '-', '_');
[data, ~] = AMPX_Naris_preprocess([],fname,'pre');
LoadExpKeys
data_remap_AMPX = AMPX_remap_sites(data, ExpKeys);


data_tsd = AMPX2tsd(data_remap_AMPX);
clear data
clear data_remap_AMPX

%% Collect the an example event

evts = all_data_pre.(fname).lg.evts;
% if exist('cfg.example', 'var') ~=1  % R061-2014_09_26 Low 
%     cfg.example = randi([1 length(evts.tstart)],1,1);
% end

ctr = mean(cat(2,evts.tstart(cfg.example),evts.tend(cfg.example)),2);
bg_tstart_idx = nearest_idx3(ctr-cfg.width/2,data_tsd.tvec);
bg_tend_idx = nearest_idx3(ctr+cfg.width/2,data_tsd.tvec);

% prepare the figure paramters
close all
h1 = figure(101);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936])

% plot the raw traces for the four corner sites
h_sub1 = subplot(3,11,[2:7, 13:18]);
color_ord = linspecer(4);
set(gca, 'ColorOrder', color_ord)
hold all
plot_loop = 1;
for ichan = [64 8 57 1]; % four corners of the probe
    %         if isempty(intersect(ichan, ExpKeys.BadChannels)) ~=1
    %             plot(trial_data.time{itrial}*1000,trial_data.trial{itrial}(ichan,:)+plot_loop*raw_plot_boost, 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
    %             text(cfg.mrk_off, median(trial_data.trial{itrial}(ichan,:))+plot_loop*raw_plot_boost, markers{plot_loop}, 'fontsize',ft_size, 'fontweight','bold','fontname', 'helvetica');
    %         else
    h.LFP = plot(data_tsd.tvec(bg_tstart_idx:bg_tend_idx),data_tsd.data(ichan, bg_tstart_idx:bg_tend_idx)+plot_loop*cfg.raw_plot_boost, 'linewidth', cfg.linewidth);
    text(data_tsd.tvec(bg_tend_idx+cfg.mrk_off), median(data_tsd.data(ichan, bg_tstart_idx:bg_tend_idx))+plot_loop*cfg.raw_plot_boost, cfg.markers{plot_loop}, 'fontsize',cfg.ft_size, 'fontweight','bold', 'fontname', 'helvetica');
    %         end
    plot_loop = plot_loop + 1;
    
end
y_lim = ylim;
% ylim([0 y_lim(2)])
xlim([data_tsd.tvec(bg_tstart_idx), data_tsd.tvec(bg_tend_idx)])
set(gca, 'xtick', data_tsd.tvec(bg_tstart_idx):.5: data_tsd.tvec(bg_tend_idx),'xTickLabel', 0:500:500, 'fontsize', cfg.ft_size, 'fontname', 'helvetica','fontweight','bold')
set(gca, 'ytick', [], 'yticklabel', [], 'TickDir','out')

rectangle('position',[evts.tstart(cfg.example), y_lim(1)+10, evts.tend(cfg.example) - evts.tstart(cfg.example), 25], 'facecolor', [.6 .6 .6], 'edgecolor', [.6 .6 .6])

% label the axes
text(data_tsd.tvec(bg_tstart_idx+round(.4*length(data_tsd.tvec(bg_tstart_idx:bg_tend_idx)))), y_lim(1)-75, 'Time (ms)','fontsize',cfg.ft_size, 'fontweight','bold', 'fontname', 'helvetica')
text(data_tsd.tvec(bg_tstart_idx-cfg.mrk_off),y_lim(1)+35, 'Ventrolateral', 'fontsize',cfg.ft_size, 'fontweight','bold', 'fontname', 'helvetica');
text(data_tsd.tvec(bg_tstart_idx-cfg.mrk_off), (plot_loop-0.45)*cfg.raw_plot_boost, 'Dorsomedial', 'fontsize',cfg.ft_size, 'fontweight','bold', 'fontname', 'helvetica');

% plot the power across the entire probe
h_sub2 = subplot(3,11,[8:11, 19:22]);
og_size = get(h_sub2, 'position');
% imagesc walkaround
all_data_pre.(fname).lg.power.power_distrib{cfg.example}(8,7)  = NaN;
nan_imagesc_ec(all_data_pre.(fname).lg.power.power_distrib{cfg.example});

% put the channel markers on the heatmap
text(1-.15,1,cfg.markers{4},'fontsize',cfg.ft_size+6, 'fontweight','bold', 'color', 'w', 'BackgroundColor', 'none', 'Margin', 1.5);
text(1-.15,8,cfg.markers{3},'fontsize',cfg.ft_size+6, 'fontweight','bold', 'color', 'w', 'BackgroundColor', 'none', 'Margin', 1.5);
text(8-.15,1,cfg.markers{2},'fontsize',cfg.ft_size+6, 'fontweight','bold', 'color', 'w', 'BackgroundColor', 'none', 'Margin', 1.5);
text(8-.15,8,cfg.markers{1},'fontsize',cfg.ft_size+6, 'fontweight','bold', 'color', 'w', 'BackgroundColor', 'none', 'Margin', 1.5);

colorbar('location','eastoutside')
set(h_sub2, 'position', og_size)
set(gca, 'xtick',[], 'ytick',[])
ax = axes('position',[0,0,1,1],'visible','off');
t_pos = get(h_sub2, 'position');
tx = text(t_pos(1),0.95,[strrep(fname, '_', '-') ' low evt: ' num2str(cfg.example)]);
    set(tx,'fontweight','bold');set(tx,'fontsize',10);
%% figure 1e: averages across rats (column) for low/high (row)
% average across the sessions for each rat.
% load('C:\temp\Naris_all_data_post.mat')
bands = {'lg', 'hg'};
for iBand = 1:2
    R061.(bands{iBand}) = mean(cat(3,all_data_pre.R061_2014_09_26.(bands{iBand}).power.power_distrib_avg, all_data_pre.R061_2014_09_27.(bands{iBand}).power.power_distrib_avg, all_data_pre.R061_2014_09_28.(bands{iBand}).power.power_distrib_avg),3);
    R054.(bands{iBand}) = mean(cat(3,all_data_pre.R054_2014_10_11.(bands{iBand}).power.power_distrib_avg, all_data_pre.R054_2014_10_12.(bands{iBand}).power.power_distrib_avg, all_data_pre.R054_2014_10_13.(bands{iBand}).power.power_distrib_avg),3);
    R049.(bands{iBand}) = mean(cat(3,all_data_pre.R049_2014_02_07.(bands{iBand}).power.power_distrib_avg, all_data_pre.R049_2014_02_08.(bands{iBand}).power.power_distrib_avg, all_data_pre.R049_2014_02_10.(bands{iBand}).power.power_distrib_avg,...
    all_data_post.R049_2014_02_07.(bands{iBand}).power.power_distrib_avg, all_data_post.R049_2014_02_08.(bands{iBand}).power.power_distrib_avg, all_data_post.R049_2014_02_10.(bands{iBand}).power.power_distrib_avg),3);
    R045.(bands{iBand}) = mean(cat(3,all_data_pre.R045_2014_04_16.(bands{iBand}).power.power_distrib_avg, all_data_pre.R045_2014_04_17.(bands{iBand}).power.power_distrib_avg, all_data_pre.R045_2014_04_18.(bands{iBand}).power.power_distrib_avg,...
        all_data_post.R045_2014_04_16.(bands{iBand}).power.power_distrib_avg, all_data_post.R045_2014_04_17.(bands{iBand}).power.power_distrib_avg, all_data_post.R045_2014_04_18.(bands{iBand}).power.power_distrib_avg),3);
end
% low gamma
h2 = figure(102);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936]);
subplot(2,4,1)
nan_imagesc_ec(R061.lg/min(min(R061.lg)));

subplot(2,4,2)
nan_imagesc_ec(R054.lg/min(min(R054.lg)));

subplot(2,4,3)
nan_imagesc_ec(R049.lg/min(min(R049.lg)));

subplot(2,4,4)
nan_imagesc_ec(R045.lg/min(min(R045.lg)));

% high gamma
set(gcf, 'position', [24 -260 1824 936]);
subplot(2,4,5)
nan_imagesc_ec(R061.hg/min(min(R061.hg)));

subplot(2,4,6)
nan_imagesc_ec(R054.hg/min(min(R054.hg)));

subplot(2,4,7)
nan_imagesc_ec(R049.hg/min(min(R049.hg)));

subplot(2,4,8)
nan_imagesc_ec(R045.hg/min(min(R045.hg)));

% make the imagescs square
axesHandles = findobj(get(h2,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
for isub = 1:8
    subplot(2,4,isub)
    set(gca, 'xtick', [], 'ytick', []);
    og_size = get(gca, 'position');
    colorbar('location', 'eastoutside');
    set(gca, 'position', og_size);
end

% add labels for each rat/gamma type
for isub = 1:4
    subplot(2,4,isub)
    t1 = title(['Rat ' num2str(isub)],  'fontsize', cfg.ft_size, 'fontweight', 'demi', 'fontname', 'helvetica');
    t_pos = get(t1, 'position');
    set(t1, 'position', [t_pos(1)-3, t_pos(2), t_pos(3)])
end

%% create the same figure with a range across all rats from the lowest to the highest

% low gamma
h3 = figure(103);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936]);
subplot(2,4,1)
nan_imagesc_ec(R061.lg);

subplot(2,4,2)
nan_imagesc_ec(R054.lg);

subplot(2,4,3)
nan_imagesc_ec(R049.lg);

subplot(2,4,4)
nan_imagesc_ec(R045.lg);

% high gamma
set(gcf, 'position', [24 -260 1824 936]);
subplot(2,4,5)
nan_imagesc_ec(R061.hg);

subplot(2,4,6)
nan_imagesc_ec(R054.hg);

subplot(2,4,7)
nan_imagesc_ec(R049.hg);

subplot(2,4,8)
nan_imagesc_ec(R045.hg);

c_max_l = max(max(cat(2,R061.lg, R054.lg, R049.lg, R045.lg)));
c_max_h = max(max(cat(2,R061.hg, R054.hg, R049.hg, R045.hg)));

c_min_l = min(min(cat(2,R061.lg, R054.lg, R049.lg, R045.lg)));
c_min_h = min(min(cat(2,R061.hg, R054.hg, R049.hg, R045.hg)));

% make the imagescs square
axesHandles = findobj(get(h3,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
for isub = 1:8
    subplot(2,4,isub)
    set(gca, 'xtick', [], 'ytick', []);
    og_size = get(gca, 'position');
    colorbar('location', 'eastoutside');
    set(gca, 'position', og_size);
    if isub <=4
    set(gca, 'clim', [c_min_l c_max_l])
    else
    set(gca, 'clim', [c_min_h c_max_h])
    end        
end

% add labels for each rat/gamma type
for isub = 1:4
    subplot(2,4,isub)
    t1 = title(['Rat ' num2str(isub)],  'fontsize', cfg.ft_size, 'fontweight', 'demi', 'fontname', 'helvetica');
    t_pos = get(t1, 'position');
    set(t1, 'position', [t_pos(1)-3, t_pos(2), t_pos(3)])
end

%% save the figures
% save the cfg.example
saveas(h1, 'D:\DATA\Paper_figs\Fig1_C_D', 'epsc')

saveas(h2, 'D:\DATA\Paper_figs\Fig1_E', 'epsc')

saveas(h3, 'D:\DATA\Paper_figs\FigS1', 'epsc')