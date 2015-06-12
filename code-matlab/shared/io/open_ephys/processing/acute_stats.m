function acute_stats(data_dir, chan_to_view, varargin)
%% Acute_stats: Take in the data specified and apply basic stats to determine the peak/trough and slope (LATER FEATURE)
%
%
%INPUTS:
%    - 'data_dir' {str}: path to file (eg:'G:\R047\2014-11-06_20-02-52')
%    - 'chan_to_view [1xN array]: list the channels to view (eg: [9 11])
%    - to use the other inputs type the name in '' then a comma, then the
%    value in the same format as it is listed in the Gather Variables
%    section.  Eg: if you want to change the decimation factor type it as:
%    ephys_data_viewer('G:\Acute\R047\HC_Stim_2014-12-04_14-15-33_r3', [9 11], 'decimate_factor', 3).
%    This can be stacked by putting commas between them Eg:
%    ephys_data_viewer('G:\Acute\R047\HC_Stim_2014-12-04_14-15-33_r3', [9 11], 'decimate_factor', 3, 'peak_jumper', 'on');

%% Gather the variables

decimate_factor = 1; % default is 3;
offset = 500;
cfg.pre_event = 0.01; % 10ms before the event. used for averaging.
cfg.post_event = .05; % 50ms post the stim 
extract_varargin
title_dir = data_dir;
cd([data_dir(1:4)])
data_dir = [data_dir '*'];
data_dir = FindFiles(data_dir);
data_dir = data_dir{1,1};
disp(['Folder to load: ' data_dir])
cd(data_dir);
listings = dir;
if isempty(FindFiles('115*'))==0 && isempty(FindFiles('110*'))
    prefix = '115_CH';
else
    prefix = '110_CH';
end


%% load the data, then decimate it
for iChan = chan_to_view
    [data.Channels{iChan, 1}.data,data.Channels{iChan, 1}.tvec,data.Channels{iChan, 1}.info] = load_open_ephys_data([prefix num2str(iChan) '.continuous']);
    data.Channels{iChan, 1}.data = decimate(data.Channels{iChan, 1}.data, decimate_factor);
    data.Channels{iChan, 1}.tvec = decimate(data.Channels{iChan, 1}.tvec, decimate_factor);
    data.Channels{iChan, 1}.info.header.sampleRate = data.Channels{iChan, 1}.info.header.sampleRate/decimate_factor;
end

disp('Loading and Decimation Complete')
%%
stat_chan = chan_to_view;
% data.Channels{stat_chan}.data(data.Channels{stat_chan}.data>1000
% [pow, pos] = findpeaks((data.Channels{stat_chan}.data), 'MINPEAKHEIGHT', 1000, 'MINPEAKDISTANCE', 10*data.Channels{iChan, 1}.info.header.sampleRate);
% 
dif = diff(data.Channels{stat_chan}.data);
dif(dif>1000) =1000;
[pow, pos] = findpeaks(dif, 'MINPEAKHEIGHT', 100, 'MINPEAKDISTANCE', 10*data.Channels{iChan, 1}.info.header.sampleRate);
% figure
% plot(dif)
% hold on
% plot(pos, pow, 'xr')
% pos2 = dif(dif>1000 & dif<1300);
% for ii = 1:length(pos2)
%         disp(ii)
%         pos3 = (find(data.Channels{stat_chan}.data == pos2(ii)));
% end
% pos = pos3;
%% Plot the channels of interest
color = linspecer(length(chan_to_view));
figHandles = get(0,'Children');
if isempty(figHandles) ==0 || sum(figHandles == 10) > 0
    close(10)
end
figure(10)
hold on
loop= 1;
for iChan = chan_to_view
    % plot(data.Channels{iChan}.tvec, data.Channels{iChan}.data+offset*(-round+1), 'Color', color(round,:))
    plot(((1:length(data.Channels{iChan}.data)).*data.Channels{iChan, 1}.info.header.sampleRate)/60, data.Channels{iChan}.data+offset*(-loop+1), 'Color', color(loop,:))
    plot(pos*data.Channels{iChan, 1}.info.header.sampleRate/60, pow, 'rx')
    leg_names{loop} =  data.Channels{iChan}.info.header.channel;
    loop = loop +1;
end
legend(leg_names)
%% extract the events
events.data = cell(length(pos),1);
events.tvec = cell(length(pos),1);
events.data_norm = cell(length(pos),1);

for ii  = 1:length(pos)
    ts = data.Channels{stat_chan}.tvec(pos(ii))-cfg.pre_event;
    [~, i_pre] =min(abs(data.Channels{stat_chan}.tvec-ts));
    tend = data.Channels{stat_chan}.tvec(pos(ii))+cfg.post_event;
    [~, i_post] =min(abs(data.Channels{stat_chan}.tvec-tend));
    events.data{ii} = data.Channels{stat_chan}.data(i_pre:i_post);
    events.tvec{ii} = data.Channels{stat_chan}.tvec(i_pre:i_post);
    events.data_norm{ii} = data.Channels{stat_chan}.data(i_pre:i_post)/(nanmean(data.Channels{stat_chan}.data(i_pre:pos(ii)-10))); 
    events.t_pre = pos(ii)-i_pre;
    events.t_post = i_post-pos(ii);
    
end
events_all = [];
for ii = 1:length(pos)
    events_all = [events_all, events.data_norm{ii}];
end
events.data_avg = nanmean(events_all,2);

%% plot all the events overlay ontop
c_ord = linspecer(length(pos));

figHandles = get(0,'Children');
if sum(figHandles == 1000) > 0
    close(1000)
end

figure(1000)
maximize
subplot(4,6,[1 2 3 7 8 9 13 14 15])
hold on
for ii = 1:length(pos)
leg{ii} = num2str(ii);
plot((1:length(events.data_norm{ii}))./mode(diff(events.tvec{ii})), events.data_norm{ii}, 'Color', c_ord(ii,:))
legend(leg)
end
% set(gca, 'ylim', [-10 10]); box on
title([title_dir '   Channel: ' chan_to_view])
subplot(4,6,[4 5 6 10 11 12 16 17 18])
 plot((events.tvec{ii}(1:events.t_pre+1000)), events.data_avg(1:events.t_pre+1000), 'k')
 hold on


%% add some stats to the plot
stats_data = events.data_avg; 
stats_data(290:325) = median(events.data_avg);
%%

% get the average peak/valley
subplot(4,6,[19:24])
box on
plot(stats_data)
peak = max(smooth(stats_data(events.t_pre:events.t_pre+1000), 5)); 
trough = min(smooth(stats_data(events.t_pre:events.t_pre+1000), 5)); 
smooth_stats_data = smooth(stats_data);
min_pow = smooth_stats_data(smooth_stats_data == trough);
min_pos = find(smooth_stats_data == trough);
[pow1, pos2] = findpeaks(smooth(stats_data(events.t_pre:events.t_pre+1000),5), 'MINPEAKHEIGHT', peak-.0001);
hold on
plot(pos2, pow1, 'rx', min_pos, min_pow-.1, 'g*')
set(gca, 'ylim', [0 2])

% subplot(6,2,11:12)
% box off; axis off
