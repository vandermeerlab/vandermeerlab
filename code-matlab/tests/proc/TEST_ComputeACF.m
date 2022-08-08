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

binsize = 10;

%% test that outputs have the expected size
cfg_acf = [];
cfg_acf.binsize = binsize;
[acf,tvec] = ComputeACF(cfg_acf, spike_data_binary); 

assert(mean(diff(tvec)) == cfg_acf.binsize); % cfg.binsive should be the same as the tvec dt
assert(tvec(ceil(end/2)) == 0) % tvec centered on 0
assert(acf(ceil(end/2)) == 0) % acf centered on 0


cfg_acf = [];
cfg_acf.binsize = binsize;
cfg_acf.maxlag = 10000; 
[acf,tvec] = ComputeACF(cfg_acf, spike_data_binary); 

assert(mean(diff(tvec)) == cfg_acf.binsize); % cfg.binsive should be the same as the tvec dt
assert(tvec(ceil(end/2)) == 0) % tvec centered on 0
assert(acf(ceil(end/2)) == 0) % acf centered on 0

cfg_acf = [];
cfg_acf.binsize = binsize;
cfg_acf.sided = 'one'; 
[acf,tvec] = ComputeACF(cfg_acf, spike_data_binary); 

assert(mean(diff(tvec)) == cfg_acf.binsize); % cfg.binsive should be the same as the tvec dt
assert(tvec(ceil(end/2)) ~= 0) % tvec not centered on 0
assert(acf(ceil(end/2)) ~= 0) % acf not centered on 0

cfg_acf = [];
cfg_acf.binsize = binsize;
cfg_acf.sided = 'onezero'; 
[acf,tvec] = ComputeACF(cfg_acf, spike_data_binary); 

assert(mean(diff(tvec)) == cfg_acf.binsize); % cfg.binsive should be the same as the tvec dt
assert(tvec(ceil(end/2)) == 0) % tvec not centered on 0
assert(acf(ceil(end/2)) == 0) % acf not centered on 0
assert(all(acf(1:ceil(end/2)) == 0)) % ensure all acf up to tvec 0 are zero


