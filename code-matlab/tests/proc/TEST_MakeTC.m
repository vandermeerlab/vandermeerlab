% setup
rng(0); % for reproducibility

% tuning variable: linear increase
tuning_dt = 1/30; 
tuning_var_tvec = 0:tuning_dt:100;
tuning_var_data = linspace(80, 400, length(tuning_var_tvec));

% spike train (logical)
spike_dt = 0.001;
spike_rate = 0.5;
spike_tvec = 0:spike_dt:100;
spike_data_binary = rand(size(spike_tvec)) > spike_rate;

% spike train (continuous)
spike_data_cont = rand(size(spike_tvec));

cfg.bins = 0:10:640;
cfg.interp = 'nearest';

%% test that outputs have the expected size
[pred_x, TC_x, x_binned] = MakeTC_1D(cfg, tuning_var_tvec, tuning_var_data, spike_tvec, spike_data_binary);
assert(length(TC_x) == length(cfg.bins) - 1); % input bins should be edges, but output corresponds to centers

[pred_x, TC_x, x_binned] = MakeTC_1D(cfg, tuning_var_tvec, tuning_var_data, spike_tvec, spike_data_cont);
assert(length(TC_x) == length(cfg.bins) - 1); % input bins should be edges, but output corresponds to centers
