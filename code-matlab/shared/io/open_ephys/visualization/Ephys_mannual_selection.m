function [events, events_to_keep, cfg] = Ephys_mannual_selection(events, cfg)
%% Ephys_mannual_selection: (Acute recording) this will bring up each event
% that has been identified using Ephys_event_detect using the data loaded
% in Ephys_load_session.  The user will be prompted to rate each event in
% order on a scale from 5:1 (one being poor).

%% define variables:



%%
c_ord = linspecer(length(cfg.events_pos));
chan_ord = linspecer(length(cfg.chan_to_view));
keep = 1;
figHandles = get(0,'Children');
if sum(figHandles == 1000) > 0
    close(1000)
end

figure(1000)
for ii = 1:length(cfg.events_pos);
    subplot(4,6,[1 2 3 7 8 9 13 14 15])
    tvec = events.tvec{ii}(:)-events.tvec{ii}(1);
    %         for ichan = cfg.chan_to_view
    plot(tvec, events.data_norm{ii}, 'Color', c_ord(ii,:))
    %         end
    vline(.025, '--g', '15ms')
    legend(num2str(ii))
    set(gca, 'ylim', [-20 20]);
    box on
    title([cfg.title_dir '   Channel: ' num2str(cfg.chan_to_view(cfg.stat_chan))])
    subplot(4, 6, [19 20])
    plot(tvec(260:400), events.data_norm{ii}(260:400), 'Color', c_ord(ii, :));
    % get user input
    user_out =input(['Keep? (y/n)\nEvent #'  num2str(ii) '\n'], 's');
    if isempty(user_out)
        user_out =input('Keep? (y/n)\n', 's');
    end
    if strcmp(user_out, 'y')
        events_to_keep.data_norm{keep} = events.data_norm{ii};
        events_to_keep.tvec{keep} = events.tvec{ii};
        events_to_keep.data{keep} = events.data{ii};
        keep = keep+1;
    end
    %     if isempty(events_to_keep)==0 ||exist('events_to_keep', 'var') ~=0
    % get the average response based on the good events only.
    if exist('events_to_keep', 'var')  && isfield(events_to_keep, 'avg') && isfield(events_to_keep, 'data_norm')
        events_to_keep.avg = nanmean(cell2mat(events_to_keep.data_norm),2);
    else
        events_to_keep.avg = events.data_norm{1};
    end
    %     end
    subplot(4,6,[4 5 6 10 11 12 16 17 18])
    tvec = events.tvec{1}(:)-events.tvec{1}(1);
    plot(tvec, events_to_keep.avg, 'k')
    set(gca, 'ylim', [-20 20]);
    vline(.025, '--g', '15ms')
    
    xlabel('Average')
    beep
end
if isfield(events_to_keep, 'tvec') ==0
    events_to_keep.tvec = tvec;
end