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
    for iChan = 1:length(cfg.chan_to_view);
        if iChan==1; h = subtightplot(6,10,[1:5 11:15]); title(cfg.title_dir)
        elseif iChan ==2; h = subtightplot(6,10,[21:25 31:35]);
        elseif iChan ==3; h = subtightplot(6,10,[6:10 16:20]);
        elseif iChan ==4; h= subtightplot(6,10,[26:30 36:40]);
        end
        tvec = events.tvec{ii}(:)-events.tvec{ii}(1);
        %         for ichan = cfg.chan_to_view
        plot(tvec, events.channel{iChan}.data_norm{ii}, 'Color', c_ord(ii,:))
        %         end
        vline([.015 .025 .035], {'--r' '--b' '--g'}, {'5ms' '15ms' '25ms'})
        legend(num2str(ii))
        set(gca, 'ylim', [-15 15], 'xlim', [0.008 .08]);
        box on
        text(.035,13,['Channel: ' num2str(cfg.chan_to_view(iChan))], 'FontSize', 16, 'Parent', h)
    end
        subtightplot(6, 10, [41:43 51:53])
        plot(tvec(260:360), events.channel{1}.data_norm{ii}(260:360), 'Color', c_ord(ii, :), 'linewidth', 4);
        
        %         tightfig
        maximize
        
        % get user input
        user_out =input(['Keep? (y/n)\nEvent #'  num2str(ii) '\n'], 's');
        if isempty(user_out)
            user_out =input('Keep? (y/n)\n', 's');
        end
        if strcmp(user_out, 'y')
            for iChan = 1:length(cfg.chan_to_view)
            events_to_keep.channel{iChan}.data_norm{keep} = events.channel{iChan}.data_norm{ii};
            end
            events_to_keep.tvec{keep} = events.tvec{ii};
            events_to_keep.data{keep} = events.data{ii};
            keep = keep+1;
        end
        %     if isempty(events_to_keep)==0 ||exist('events_to_keep', 'var') ~=0
        % get the average response based on the good events only.
        if exist('events_to_keep', 'var')  && isfield(events_to_keep.channel{1,1}, 'avg') && isfield(events_to_keep.channel{1,1}, 'data_norm')
            for iChan = 1:length(cfg.chan_to_view)
                events_to_keep.channel{iChan}.avg = nanmean(cell2mat(events_to_keep.channel{iChan}.data_norm),2);
            end
        elseif exist('events_to_keep', 'var') 
            for iChan = 1:length(cfg.chan_to_view)
                events_to_keep.channel{iChan}.avg = events.channel{iChan}.data_norm{1};
            end
        end
        %     end
        if exist('events_to_keep', 'var')  && isfield(events_to_keep.channel{1,1}, 'avg')
            
            subtightplot(6,10, [44:50 54:60])
            tvec = events.tvec{1}(:)-events.tvec{1}(1);
            plot(tvec, events_to_keep.channel{iChan}.avg, 'k')
            set(gca, 'ylim', [-15 15], 'xlim', [.008 .04]);
            vline([.015 .025 .035], {'--r' '--b' '--g'}, {'5ms' '15ms' '25ms'})
            
            xlabel('Average')
        end
end
if isfield(events_to_keep, 'tvec') ==0
    events_to_keep.tvec = tvec;
end