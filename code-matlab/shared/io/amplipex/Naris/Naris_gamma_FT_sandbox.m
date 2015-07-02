%% Naris to FT trail sandbox

%% load the events file and get ivs

% defaults
gamma_evt_window = 0.05;

%% load data

cfg.fname = 'R060-2015-01-31';
exp = 'pre';
cfg.df = 10;
cfg.low_gamma= [40 55];
cfg.high_gamma = [70 85];
[cfg] = Naris_cfgs(cfg);
if strcmp(cfg.fname(1:4), 'R053') || strcmp(cfg.fname(1:4), 'R060')
    cfg.tts_to_process = [5:13,15];
    cfg.Naris_exp = {'pre', 'right', 'left', 'post'};
    tetrodes = tts(cfg);
    cfg.tetrodes = tetrodes.channels_array(cfg.chan);
    cfg.channels = tetrodes.channels_array;
end
data = AMPX_loadData([cfg.fname '-' exp '.dat'], cfg.channels(1:4),cfg.df); % this is where you can specify which tt to load. 

%% convert to FT
[data_ft] = AMPX_makeft(data);

%% detect gamma events

[low_gamma_iv, high_gamma_iv, random_lg_iv, random_hg_iv, ~, ~] = AMPX_Gamma_detect_nsb2014 (data, 'channel_of_interest', 1, 'display', 'off');

%% Split data_FT into gamma events
% restructure ivs and correct for Fs
low_events = ([low_gamma_iv.tstart, low_gamma_iv.tend])*data_ft.hdr.Fs;
high_events = ([high_gamma_iv.tstart, high_gamma_iv.tend])*data_ft.hdr.Fs;

%trial split into FT format
    [data_trl_low, ~, ~] = AMPX_trial_split(data_ft, low_events, gamma_evt_window); % this may not be appropriote since it creates trials of the same length. 
    [data_trl_high, ~, ~] = AMPX_trial_split(data_ft, high_events, gamma_evt_window);

