function [cfg] = Ephys_stats(events_to_keep, cfg)
% Ephys_stats: gathers basic statistics for each event following the
% loading, event detection and event selection/inspection phases on the
% Ephys pipeline

%% define the defaults
figHandles = get(0,'Children');
if sum(figHandles == 2000) > 0
    close(2000)
end
figure(2000)
tvec = events_to_keep.tvec{1}(:)-events_to_keep.tvec{1}(1);
hold on
plot(tvec, events_to_keep.avg, 'k')
    set(gca, 'ylim', [-10 10]);
    vline([0.015 .025 0.065 0.075], {'--r', '--g', '--r', '--g'}, {'5ms', '15ms', 'P2 5ms', 'P2 15ms'})
    
    title([num2str(length(events_to_keep.data_norm)) '/' num2str(length(cfg.events_pos))])
%% get the inputs
% [x, y] = ginput(2);
% 
% x_1 = find(round(tvec * 10000) / 10000 ==round(x(1) * 10000) / 10000);
% x_2 = find(round(tvec * 10000) / 10000 ==round(x(2) * 10000) / 10000);
% range = x_1(1) : x_2(1);
% 
% data_range = events_to_keep.avg(range);
% 
% [stats.max.val, stats.max_idx ]  = max(smooth(data_range, 5)); 

end
