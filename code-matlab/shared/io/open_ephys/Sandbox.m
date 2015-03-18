%% acute inspector sandbox

%% load and process some data
data_dir = 'G:\Acute\';
fname = 'R059-2';
cfg.decimate_factor = 1;
cfg.chan_to_view = [9 10];
[data, cfg] = Ephys_load_session(data_dir, fname, cfg);

%% plot all the data channels
cfg.offset = 500;

Ephys_chan_plot(data, cfg)


%% run the event detection
cfg.stat_chan = 9;
[data, events, cfg] = Ephys_event_detect(data, cfg);

%% replot the data with the events

Ephys_chan_plot(data, cfg)

%% run the mannual selection script

[~,events_to_keep, cfg]= Ephys_mannual_selection(events, cfg);

%% get some stats (under development)


