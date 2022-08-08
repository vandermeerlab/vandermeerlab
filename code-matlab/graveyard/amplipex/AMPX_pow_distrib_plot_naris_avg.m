function AMPX_pow_distrib_plot_naris_avg(power, trial_data, data, varargin)
%AMPX_pow_distrib_plot_naris
%  This is used to make specific figures for the Naris paper.
%  it only uses 4 corner channels and cleans things up a bit.  nothing
%  under the hood has changed.
%
%
%  Plot the power distributions from AMPX_pow_distrib
%
%  This will generate plots of the raw signal for the windowsize used to
%  compute the power distribution with the filtered signal over top
%  (optional), the phase of the signal and the imagesc of the power in the
%  corrrect probe layout.
%
%  Saving is optional.
% EC v2 2015-01-07
%% extract the input parameters and create some using the ExpKeys
run(FindFile('*Keys.m'));
save_fig = 'on';
min_freq = power.params.freq(1);
max_freq = power.params.freq(2);
if strcmp(data.hdr.Filename(end-7:end-4), 'post')
    session_type = 'post';
elseif strcmp(data.hdr.Filename(end-6:end-4), 'pre')
    session_type = 'pre';
else
    session_type = 'task';
end
% butter_order = 10;
raw_plot_boost = 150; % use a rescale function to see the data in a different way.  use this as a flag
if strcmp(ExpKeys.ProbeType, 'A8x8')
    psd_sub_order_r = [37 45 53 61 5 13 21 29 36 44 52 60 4 12 20 28 38 46 54 62 6 14 22 30 35 43 51 59 3 11 19 27 39 47 55 63 7 15 23 31 34 42 50 58 2 10 18 26 40 48 56 64 8 16 24 32 33 41 49 57 1 9 17 25];
    %     mid_chan = [35 25 16 52 39 30 6 64]; %
    %     mid_chan = [64  60 47 43 ];
    mid_chan = [57 61 34 38 ];
    %     mid_chan = [ 57 2 22 43];
elseif strcmp(ExpKeys.ProbeType, 'Buz64')
    psd_sub_order_r = [33 41 49 57 1 9 17 25 40 48 56 64 8 16 24 32 34 42 50 58 2 10 18 26 39 47 55 63 7 15 23 31 35 43 51 59 3 11 19 27 38 46 54 62 6 14 22 30 36 44 52 60 4 12 20 28 37 45 53 61 5 13 21 29];
    mid_chan = [48 64 15 31 44 60 12 28];
end
ft_size = 20;
extract_varargin

if min_freq<65 && min_freq> 20; type = 'low_gamma'; clim = [0 900]; type_title = 'low gamma'; raw_plot_boost = 200;line_end  =.1; % use a rescale function to see the data in a different way.  use this as a flag
elseif min_freq>66 && max_freq<100; type = 'high_gamma'; type_title = 'high gamma'; clim = [10 300]; raw_plot_boost = 150; line_end  =.1;
elseif min_freq>5 && max_freq <10; type = 'theta'; type_title = 'theta';  clim = [50 500]; raw_plot_boost = 400; line_end  =1;
elseif min_freq>5 && max_freq <15; type = 'spindles'; type_title = 'spindles';  clim = [50 500]; raw_plot_boost = 400; line_end  =1;
elseif 100<=min_freq; type = 'high_freq'; type_title = 'high freq with randomized iei'; clim = [0.5 2.5];
elseif min_freq>1 && max_freq <5; type = 'delta'; type_title = 'delta';  clim = [50 500]; raw_plot_boost = 400; line_end  =1;
    
end
% if strcmp(ExpKeys.Subject, 'R061'); raw_plot_boost = raw_plot_boost/3; end
line_end = line_end*1000;
if strcmp(ExpKeys.Subject, 'R045')==1
    %     mid_chan = [38 27 14 61 42 19 1 58];
end
%% Get the filtered data.
if strcmp(type, 'high_gamma')
    mrk_off = -5;
else
    mrk_off = -50;
end
% [data_ft_filt, data_filtered] = AMPX_filter(data, min_freq-10 , max_freq+10, butter_order);
% [data_trl_filt, ~, ~] = AMPX_trial_split(data_ft_filt, trial_data.cfg.trl, power.params.window_size);

%% create the power figures using the raw data,
    %%
    figure(1000);
    maximize
    subplot(1,6,1:2)
    %     color_ord = get(gca,'colororder');
    color_ord = linspecer(4);
    %     color_ord = repmat([0 0 0],4,1);
    %     color_ord =  color_ord(1:3,:);
    %     color_ord = [color_ord; [0 .75 .75]];
    set(gca, 'ColorOrder', color_ord)
    hold all
    plot_loop = 1;
%     for ichan = mid_chan
        markers = {'#', '+', 'x', 'o'};
%         %         chan_ind = find(ichan == ExpKeys.Probe_layout);
%         if isempty(intersect(ichan, ExpKeys.BadChannels)) ~=1
%             plot(trial_data.time{itrial}*1000,trial_data.trial{itrial}(ichan,:)+plot_loop*raw_plot_boost, 'Color', [0.5 0.5 0.5], 'LineStyle', '--')
%             text(mrk_off, median(trial_data.trial{itrial}(ichan,:))+plot_loop*raw_plot_boost, markers{plot_loop}, 'fontsize',ft_size, 'fontweight','bold','fontname', 'helvetica');
%             
%         else
%             plot(trial_data.time{itrial}*1000,trial_data.trial{itrial}(ichan,:)+plot_loop*raw_plot_boost)
%             text(mrk_off, median(trial_data.trial{itrial}(ichan,:))+plot_loop*raw_plot_boost, markers{plot_loop}, 'fontsize',ft_size, 'fontweight','bold', 'fontname', 'helvetica');
%         end
%         %         plot(data_trl_filt.time{itrial}, data_trl_filt.trial{itrial}(ichan,:)+ plot_loop*raw_plot_boost, ':k')
%         plot_loop = plot_loop + 1;
%         %         if plot_loop == 3
%         %             line(0:line_end/100:line_end, plot_loop*raw_plot_boost)
%         %             plot_loop = plot_loop +1;
%         %         end
%     end
    %%
    %     legend ('64', '15', '23', '47', '60', '5', '20', '43', 'Location', 'NorthWestOutside')% '45', '20', '5', '60', '47', '23', '15', '64')
    xlim([0 trial_data.time{1}(end)*1000]); ylim([-0.5*raw_plot_boost (plot_loop+.05)*raw_plot_boost]);
    set(gca,  'xtick',[0 line_end] , 'ytick', []);
    set (gcf, 'color', [1 1 1])
    color = get(gca,'Color');
    set(gca,'XColor',color,'YColor',color,'TickDir','out')
%     %     xlabel('Time (ms)', 'FontSize', 16, 'FontWeight','bold')
%     hline(0, 'k')
%     if strcmp(type, 'high_gamma')
%         text(-0, -15, num2str(0),'fontsize',ft_size-4, 'fontweight','bold', 'fontname', 'helvetica');
%         text((trial_data.time{1}(end)*1000)-5, -15, num2str(trial_data.time{1}(end)*1000),'fontsize',ft_size-4, 'fontweight','bold', 'fontname', 'helvetica');
%         text(trial_data.time{1}(end)*1000-0.5*trial_data.time{1}(end)*1000-5, -25, 'Time (ms)','fontsize',ft_size-4, 'fontweight','bold', 'fontname', 'helvetica')
%         text(.035, 25, 'Ventrolateral', 'fontsize',ft_size, 'fontweight','bold', 'fontname', 'helvetica');
%         %     text(.035, (plot_loop+.5)*raw_plot_boost/2.5, 'Dorsal Lateral', 'fontsize',16, 'fontweight','bold');
%         %     text(.035, (plot_loop+.5)*raw_plot_boost/2, 'Ventral Medial', 'fontsize',16, 'fontweight','bold');
%         text(.035, (plot_loop-0.3)*raw_plot_boost, 'Dorsomedial', 'fontsize',ft_size, 'fontweight','bold', 'fontname', 'helvetica');
%         
%     else
%         text(-10, -25, num2str(0),'fontsize',ft_size-4, 'fontweight','bold', 'fontname', 'helvetica');
%         text((trial_data.time{1}(end)*1000)-50, -25, num2str(trial_data.time{1}(end)*1000),'fontsize',ft_size-4, 'fontweight','bold', 'fontname', 'helvetica');
%         text(trial_data.time{1}(end)*1000-0.5*trial_data.time{1}(end)*1000-100, -65, 'Time (ms)','fontsize',ft_size-4, 'fontweight','bold', 'fontname', 'helvetica')
%         text(.035, 45, 'Ventrolateral', 'fontsize',ft_size, 'fontweight','bold', 'fontname', 'helvetica');
%         %     text(.035, (plot_loop+.5)*raw_plot_boost/2.5, 'Dorsal Lateral', 'fontsize',16, 'fontweight','bold');
%         %     text(.035, (plot_loop+.5)*raw_plot_boost/2, 'Ventral Medial', 'fontsize',16, 'fontweight','bold');
%         text(.035, (plot_loop-0.3)*raw_plot_boost, 'Dorsomedial', 'fontsize',ft_size, 'fontweight','bold', 'fontname', 'helvetica');
%     end
    %%
    subplot(1,5,3:5)
    %             set(gcf, 'Renderer', 'opengl')
    %         x = imagesc(power.power_distrib{itrial});
    %         set(x, 'AlphaData', ~isnan(power.power_distrib{itrial}))
    %     set(gcf, 'Renderer', 'painters')
    % imagesc walkaround
    h = nan_imagesc_ec(power.power_distrib_avg);
    %
    plt_ft_size = 42;
%     text(1-.15,1,markers{4},'fontsize',plt_ft_size, 'fontweight','bold', 'color', 'k', 'BackgroundColor', 'none', 'Margin', 1.5);
%     text(1-.15,8,markers{3},'fontsize',plt_ft_size, 'fontweight','bold', 'color', 'k', 'BackgroundColor', 'none', 'Margin', 1.5);
%     text(8-.15,1,markers{2},'fontsize',plt_ft_size, 'fontweight','bold', 'color', 'k', 'BackgroundColor', 'none', 'Margin', 1.5);
%     text(8-.15,8,markers{1},'fontsize',plt_ft_size, 'fontweight','bold', 'color', 'k', 'BackgroundColor', 'none', 'Margin', 1.5);
    %     text(2-.2,7,num2str(mid_chan(5)),'fontsize',16, 'fontweight','bold',  'color', color_ord(1,:), 'BackgroundColor', [1 1 1], 'Margin', 1);
    %     text(4-.2,7,num2str(mid_chan(6)),'fontsize',16, 'fontweight','bold',  'color', color_ord(2,:), 'BackgroundColor', [1 1 1], 'Margin', 1);
    %     text(6-.2,7,num2str(mid_chan(7)),'fontsize',16, 'fontweight','bold',  'color', color_ord(3,:), 'BackgroundColor', [1 1 1], 'Margin', 1);
    %     text(8-.2,7,num2str(mid_chan(8)),'fontsize',16, 'fontweight','bold',  'color', color_ord(4,:), 'BackgroundColor', [1 1 1], 'Margin', 1);
    colorbar('location','southoutside')
    %     xlabel([ type_title  ' ' num2str(min_freq) 'Hz to ' num2str(max_freq) 'Hz'], 'fontsize', ft_size, 'fontweight','bold', 'fontname', 'helvetica')
    set(gca, 'xtick',[], 'ytick',[])
    ax = axes('position',[0,0,1,1],'visible','off');
    tx = text(0.46,0.95,[data.hdr.Filename(5:end-4) ' ' type_title ' avg']);
    set(tx,'fontweight','bold');set(tx,'fontsize',10);
    %     title([num2str(min_freq) 'Hz to ' num2str(max_freq) 'Hz'])
    set(gcf,'Units','pixels','Position',[0 0 1920 1080])
    %%
    %     if strcmp(save_fig,'on'); print(gcf,'-dpng','-r300',['D:\DATA\' data.hdr.Filename(5:8) '\' data.hdr.Filename(5:19) '\' date '\' session_type '\' type '_pow_figures\' data.hdr.Filename(5:end-4) '_' num2str(itrial) '_spec.png']);end
    mkdir(['D:\DATA\' data.hdr.Filename(5:8) '\' data.hdr.Filename(5:19) '\' date '\' session_type '\' type '_pow_figures\'])
    target_dir = ['D:\DATA\' data.hdr.Filename(5:8) '\' data.hdr.Filename(5:19) '\' date '\' session_type '\' type '_pow_figures\' data.hdr.Filename(5:end-4) '_avg_spec.png'];
    export_fig(target_dir, '-native')
    if strcmp(save_fig,'on'); saveas(gcf,['D:\DATA\' data.hdr.Filename(5:8) '\' data.hdr.Filename(5:19) '\' date '\' session_type '\' type '_pow_figures\' data.hdr.Filename(5:end-4) '_avg_spec.fig']);end
    naris_dir =     ['G:\Naris\paper_figs\' data.hdr.Filename(5:end-4) '_avg_spec.png'];
    export_fig(naris_dir, '-native')
    close(1000)

end

