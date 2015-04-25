function [data, events, cfg] = Ephys_event_detect(data, cfg)
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

cfg.pre_event = 0.01; % 10ms before the event. used for averaging.
cfg.post_event = .08; % 50ms post the stim
if isfield(cfg, 'stat_chan') == 0
    warning('You need to specify  channel to use for event detection as cfg.stat_chan ...setting to default first channel')
    cfg.stat_chan = 1;
else
    cfg.stat_chan = find(data.labels == cfg.stat_chan);
end

dif = diff(data.Channels{cfg.stat_chan}.data);
dif(dif>1000) =1000;
[pow, pos] = findpeaks(dif, 'MINPEAKHEIGHT', 100, 'MINPEAKDISTANCE', 30000);%.04*data.Channels{cfg.stat_chan, 1}.info.header.sampleRate);
% pos = pos(1:2:end);
% pow = pow(1:2:end);
%% extract the events
events.data = cell(length(pos),1);
events.tvec = cell(length(pos),1);
events.channel = cell(length(cfg.chan_to_view),1);
for iChan = 1:length(cfg.chan_to_view)
    events.channel{iChan}.data_norm = cell(length(pos),1);
end

for ii  = 1:length(pos)
    ts = data.Channels{cfg.stat_chan}.tvec(pos(ii))-cfg.pre_event;
    [~, i_pre] =min(abs(data.Channels{cfg.stat_chan}.tvec-ts));
    tend = data.Channels{cfg.stat_chan}.tvec(pos(ii))+cfg.post_event;
    [~, i_post] =min(abs(data.Channels{cfg.stat_chan}.tvec-tend));
    events.data{ii} = data.Channels{cfg.stat_chan}.data(i_pre:i_post);
    events.tvec{ii} = data.Channels{cfg.stat_chan}.tvec(i_pre:i_post);
    for iChan = 1:length(cfg.chan_to_view);
        events.channel{iChan}.data_norm{ii} = data.Channels{iChan}.data(i_pre:i_post)/(nanmean(data.Channels{iChan}.data(i_pre:pos(ii)-10)));
    end
    events.t_pre = pos(ii)-i_pre;
    events.t_post = i_post-pos(ii);
    
end
events_all = [];
for ii = 1:length(pos)
    for iChan = 1:length(cfg.chan_to_view);
    events.all = [events_all, events.channel{iChan}.data_norm{ii}];
    end
end
events.data_avg = nanmean(events_all,2);
cfg.events_pos = pos;
cfg.events_pow = pow;