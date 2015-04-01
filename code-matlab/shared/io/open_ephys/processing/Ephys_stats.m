function [cfg] = Ephys_stats(events_to_keep, cfg)
% Ephys_stats: gathers basic statistics for each event following the
% loading, event detection and event selection/inspection phases on the
% Ephys pipeline

%% param,
if isfield(cfg, 'x_lim')~=1
    cfg.x_lim = .035;
end
    

%% define the defaults
figHandles = get(0,'Children');
if sum(figHandles == 2000) > 0
    close(2000)
end
figure(2000)

maximize
tvec = events_to_keep.tvec{1}(:)-events_to_keep.tvec{1}(1);
hold on
for iChan = 1:length(cfg.chan_to_view)
        if iChan==1; h = subtightplot(2,2,1); 
        elseif iChan ==2; h = subtightplot(2,2,3);
        elseif iChan ==3; h = subtightplot(2,2,2);
        elseif iChan ==4; h= subtightplot(2,2,4);
        end
        
plot(tvec, events_to_keep.channel{iChan}.avg, 'k')
    set(gca, 'ylim', [-7 8], 'xlim', [0.0 cfg.x_lim]);
    vline([0.015 .025 0.065 0.075], {'--r', '--g', '--r', '--g'}, {'5ms', '15ms', 'P2 5ms', 'P2 15ms'})
    text(0.0025,2.5,['Channel: ' num2str(cfg.chan_to_view(iChan))], 'FontSize', 16, 'Parent', h)
    box off
    set(gca,'YTickLabel',[],'XTickLabel',[])
end
text(0.0025, 7,[num2str(length(events_to_keep.channel{iChan}.data_norm)) '/' num2str(length(cfg.events_pos))], 'FontSize', 22)

%%
figure(3000)

maximize
tvec = events_to_keep.tvec{1}(:)-events_to_keep.tvec{1}(1);
hold on
for iChan = 1:length(cfg.chan_to_view)
        if iChan==1; h = subtightplot(2,2,1); 
        elseif iChan ==2; h = subtightplot(2,2,3);
        elseif iChan ==3; h = subtightplot(2,2,2);
        elseif iChan ==4; h= subtightplot(2,2,4);
        end
        
plot(tvec, events_to_keep.channel{iChan}.avg, 'k')
    set(gca, 'ylim', [-7 8], 'xlim', [0.005 .015]);
    vline([0.015 .025 0.065 0.075], {'--r', '--g', '--r', '--g'}, {'5ms', '15ms', 'P2 5ms', 'P2 15ms'})
    text(0.0025,2.5,['Channel: ' num2str(cfg.chan_to_view(iChan))], 'FontSize', 16, 'Parent', h)
    box off
    set(gca,'YTickLabel',[],'XTickLabel',[])
end
text(0.01, 7,[num2str(length(events_to_keep.channel{iChan}.data_norm)) '/' num2str(length(cfg.events_pos))], 'FontSize', 22)
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
