function AMPX_Naris_fig_2_task(all_data_pre, all_data_post, varargin)
%% AMPX_Naris_fig_2_task: generates figure 2 which contains the heatmaps
% for the average rewarded/approach gamma 50 and gamma 80 events for the
% two rats that reached task criterion
%
%
% EC - 2016-05-31

%% extract inputs
extract_varargin

%% load the data for the task events
load('C:\temp\Naris_all_data_task')

%% R049
h1 = figure(21);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936]);

Session_list = fieldnames(all_data_task);
bands = {'lg_reward', 'hg_reward', 'lg_approach', 'hg_approach'};
for iSess = 1:length(Session_list)
    for  iBand = 1:length(bands)
        if isfield(all_data_task.(Session_list{iSess}).(bands{iBand}).power, 'power_avg') || isempty(all_data_task.(Session_list{iSess}).(bands{iBand}).power.power_distrib_avg)
            all_data_task.(Session_list{iSess}).(bands{iBand}).power.power_distrib_avg = NaN*zeros(8,8);
            warning(['Session ' Session_list{iSess} '  ' bands{iBand} '  did not contain any events'])
        end
    end
end
subplot(2,4,1)
R049.low.approach = nanmean(cat(3,all_data_task.R049_2014_02_07.lg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_08.lg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_10.lg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.low.approach/min(min(R049.low.approach)))

subplot(2,4,2)
R049.low.reward = nanmean(cat(3,all_data_task.R049_2014_02_07.lg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_08.lg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_10.lg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.low.reward/min(min(R049.low.reward)))

subplot(2,4,5)
R049.high.approach = nanmean(cat(3,all_data_task.R049_2014_02_07.hg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_08.hg_approach.power.power_distrib_avg, all_data_task.R049_2014_02_10.hg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.high.approach/min(min(R049.high.approach)))

subplot(2,4,6)
R049.high.reward = nanmean(cat(3,all_data_task.R049_2014_02_07.hg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_08.hg_reward.power.power_distrib_avg, all_data_task.R049_2014_02_10.hg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R049.high.reward/min(min(R049.high.reward)))

% R045  (R4)
figure(21)
subplot(2,4,3)
R045.low.approach = nanmean(cat(3,all_data_task.R045_2014_04_16.lg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_17.lg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_18.lg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.low.approach/min(min(R045.low.approach)))

subplot(2,4,4)
R045.low.reward = nanmean(cat(3,all_data_task.R045_2014_04_16.lg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_17.lg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_18.lg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.low.reward/min(min(R045.low.reward)))

subplot(2,4,7)
R045.high.approach = nanmean(cat(3,all_data_task.R045_2014_04_16.hg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_17.hg_approach.power.power_distrib_avg, all_data_task.R045_2014_04_18.hg_approach.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.high.approach/min(min(R045.high.approach)))

subplot(2,4,8)
R045.high.reward = nanmean(cat(3,all_data_task.R045_2014_04_16.hg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_17.hg_reward.power.power_distrib_avg, all_data_task.R045_2014_04_18.hg_reward.power.power_distrib_avg), 3);
nan_imagesc_ec(R045.high.reward/min(min(R045.high.reward)))

axesHandles = findobj(get(h1,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
for isub = 1:8
    subplot(2,4,isub)
    set(gca, 'xtick', [], 'ytick', []);
    og_size = get(gca, 'position');
    colorbar('location', 'eastoutside');
    set(gca, 'position', og_size);
end

%% export part 1
if exist('save_fig')
    saveas(h1, 'D:\DATA\Paper_figs\Fig2_A', 'epsc')
    saveas(h1, 'D:\DATA\Paper_figs\Fig2_A', 'fig')
end
%% get the counts across all pre and post recording sessions
clear all_data_task
% 
% load('C:\temp\Naris_all_data_pre.mat')
% load('C:\temp\Naris_all_data_post.mat')

%% get R^2 across sessions
sess_list = fieldnames(all_data_pre);
bands = {'lg','lg_ran', 'hg', 'hg_ran'};% 'hg'};
all.lg.rsq = []; all.hg.rsq = [];all.lg_ran.rsq = []; all.hg_ran.rsq = [];
all.lg.rsq2 = []; all.hg.rsq2 = [];all.lg_ran.rsq2 = []; all.hg_ran.rsq2 = [];

for iSess = 1:length(sess_list)
    for iband = 1:length(bands)
        all.(bands{iband}).rsq = [all.(bands{iband}).rsq, all_data_pre.(sess_list{iSess}).(bands{iband}).power.stats.rsq];
        all.(bands{iband}).rsq2 = [all.(bands{iband}).rsq2, all_data_pre.(sess_list{iSess}).(bands{iband}).power.stats.rsq2];
        if strcmp(sess_list{iSess}(1:4), 'R045') %|| strcmp(sess_list{iSess}(1:4), 'R049')
            all.(bands{iband}).rsq = [all.(bands{iband}).rsq, all_data_post.(sess_list{iSess}).(bands{iband}).power.stats.rsq];
            all.(bands{iband}).rsq2 = [all.(bands{iband}).rsq2, all_data_post.(sess_list{iSess}).(bands{iband}).power.stats.rsq2];
        end
        
    end
end
%%
[~, plane_stats.low] = AMPX_plane_count_hist(all.lg, all.lg_ran, 'low');
hax1 = gca;
hf1 = figure(22);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'position', [24 -260 1824 936]);
s1 = subplot(5,2,[1 3]);
pos=get(s1,'Position');
delete(s1);
hax2=copyobj(hax1,hf1);
set(hax2, 'Position', pos);


[~, plane_stats.high] = AMPX_plane_count_hist(all.hg, all.hg_ran, 'high');
hax11 = gca;
hf1 = figure(22);
s2 = subplot(5,2,[2 4]);
pos=get(s2,'Position');
delete(s2);
hax2=copyobj(hax11,hf1);
set(hax2, 'Position', pos);

axesHandles = findobj(get(hf1,'Children'), 'flat','Type','axes');
set(axesHandles, 'xtick', 0:20:80, 'ytick', 0:50:500, 'xlim', [0 90], 'ylim', [0 350]);
og_size = get(axesHandles, 'position');
set(axesHandles, 'position', og_size);
set(axesHandles, 'fontsize', 18, 'fontname', 'helvetica')


%% export part 1
if exist('save_fig')
    saveas(hf1, 'D:\DATA\Paper_figs\Fig2_B_C', 'epsc')
    saveas(hf1, 'D:\DATA\Paper_figs\Fig2_B_C', 'fig')
end
