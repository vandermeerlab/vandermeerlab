% setup -- this section SHOULD rsult in shared variables that all tests can
% access, but for some reason doesn't work in R2018a
%
%empty_ts = ts();
%
%one_ts_data = (0:1:10)';
%one_ts = ts; one_ts.t{1} = one_ts_data; one_ts.label{1} = 'label1';
%
%rvec = [0 10]; rvec_iv = iv(rvec);

%% test if format is preserved with raw input times
empty_ts = ts();
rvec = [0 10];

empty_ts_out = restrict(empty_ts,rvec(1),rvec(2));
assert(CheckTS(empty_ts_out) == 1);

%% test if format is preserved with iv input times
empty_ts = ts();
rvec = [0 10]; rvec_iv = iv(rvec);

empty_ts_out = restrict(empty_ts,rvec_iv);
assert(CheckTS(empty_ts_out) == 1);

%% test if ts functionality works as expected
one_ts_data = (0:1:10)';
one_ts = ts; one_ts.t{1} = one_ts_data; one_ts.label{1} = 'label1';

rvec_1 = [2.5, 8.5]; rvec_1_iv = iv(rvec_1);

one_ts_out = restrict(one_ts, rvec_1(1), rvec_1(2));
assert(all(one_ts_out.t{1} == (3:1:8)') == 1);

one_ts_out = restrict(one_ts, rvec_1_iv);
assert(all(one_ts_out.t{1} == (3:1:8)') == 1);
