% setup
rng(0); % for reproducibility

% Sanity Check: Check the if the output tvec matches the input specifcations
input_ts = sort(rand(1,100));
max_lag = 1; % in seconds
binsize = 0.001; % in seconds
cfg_xc.max_t = max_lag;%Maximum lag at which cross-correlation must be calculated.
cfg_xc.binsize = binsize;
[output_xc, output_tvec] = ccf(cfg_xc, input_ts, input_ts);
% Max value of output_tvec should be max_t
assert(max(output_tvec) == max_lag);
% See if the output_tvec is of the correct length
assert(length(output_tvec) == length(-1*max_lag:binsize:1*max_lag));

% Test if the peak in ccf occurs at lag 0 if the same timeseries is passed
% twice (while calculating autocorrelation

input_ts = sort(rand(1,1000));
max_lag = 0.5; % in seconds
binsize = 0.001;
cfg_xc.max_t = max_lag;
cfg_xc.binsize = binsize;
[output_xc, output_tvec] = ccf(cfg_xc, input_ts, input_ts);
[~, max_idx] = max(output_xc);
assert(output_tvec(max_idx) == 0);

% Test of the peak in ccf occurs at the expected delay when the two
% timestamps are at a fixed delay
fixed_delay = 0.25;
input_ts1 = sort(rand(1,1000));
input_ts2 = input_ts1 + fixed_delay;
max_lag = 0.5; % in seconds
binsize = 0.001;
cfg_xc.binsize = binsize;
[output_xc, output_tvec] = ccf(cfg_xc, input_ts1, input_ts2);
[~, max_idx] = max(output_xc);
assert(output_tvec(max_idx) == fixed_delay || ...
    output_tvec(max_idx) == -1*fixed_delay);

% Test if a spike train with theta modulated firing has max value at time
% lags corresponding to theta cycle length, not inculding the lag at 0
input_ts = 0:0.125:10;
max_lag = 1; % in seconds
binsize = 0.001;
cfg_xc.max_t = max_lag;
cfg_xc.binsize = binsize;
[output_xc, output_tvec] = ccf(cfg_xc, input_ts, input_ts);
right_half = output_xc(output_tvec > 0);
right_tvec = output_tvec(output_tvec > 0);
[~, max_idx] = max(right_half);
assert(right_tvec(max_idx) == 0.125);


