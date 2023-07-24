function [all_data, cfg] = Ephys_acute_preprocess(fname, data_dir)

%% load and process some data
% data_dir = 'G:\Acute\';
% fname = 'R058-66';
cfg.decimate_factor = 1;
cfg.chan_to_view = [1 5 25 29];
[data, cfg] = Ephys_load_session(data_dir, fname, cfg);

% %% plot all the data channels
% cfg.offset = 500;
% 
% Ephys_chan_plot(data, cfg)


% run the event detection
cfg.stat_chan = 5;
[data, events, cfg] = Ephys_event_detect(data, cfg);

%% replot the data with the events

Ephys_chan_plot(data, cfg)

%% run the mannual selection script

[~,events_to_keep, cfg]= Ephys_mannual_selection(events, cfg);

%% get some stats (under development)

[cfg] = Ephys_stats(events_to_keep, cfg);

%% save both figures and the strcut
all_data.events = events;
all_data.events_to_keep = events_to_keep;
all_data.data = data;
all_data.cfg = cfg;

save(['G:\Acute\' fname(1:4) '\all_data\' strrep(fname, '-', '_') '_all_data'], 'all_data', '-v7.3')

print(2000, ['G:\Acute\' fname(1:4) '\all_data\' strrep(fname, '-', '_') '_all_data.png']);
saveas(2000, ['G:\Acute\' fname(1:4) '\all_data\' strrep(fname, '-', '_') '_all_data.fig']);
print(3000, ['G:\Acute\' fname(1:4) '\all_data\' strrep(fname, '-', '_') '_all_data_closeup.png']);
saveas(3000, ['G:\Acute\' fname(1:4) '\all_data\' strrep(fname, '-', '_') '_all_data_closeup.fig']);

close all
end