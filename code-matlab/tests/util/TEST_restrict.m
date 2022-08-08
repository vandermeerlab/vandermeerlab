% setup -- this section SHOULD rsult in shared variables that all tests can
% access, but for some reason doesn't work in R2018a
%
%empty_ts = ts();
%
%one_ts_data = (0:1:10)';
%one_ts = ts; one_ts.t{1} = one_ts_data; one_ts.label{1} = 'label1';
%
%rvec = [0 10]; rvec_iv = iv(rvec);

%% preserve ts data type with raw input times
empty_ts = ts();
rvec = [0 10];

% start-end input argument
empty_ts_out = restrict(empty_ts, rvec(1), rvec(2));
assert(CheckTS(empty_ts_out) == 1);

% using iv input argument
rvec_iv = iv(rvec);
empty_ts_out = restrict(empty_ts, rvec_iv);
assert(CheckTS(empty_ts_out) == 1);

%% ts functionality
one_ts_data = (0:1:10)';
one_ts = ts; one_ts.t{1} = one_ts_data; one_ts.label{1} = 'label1';

rvec_1 = [2.5, 8.5]; rvec_1_iv = iv(rvec_1);

% basic, non-edge case (start-end input arguments)
one_ts_out = restrict(one_ts, rvec_1(1), rvec_1(2));
assert(all(one_ts_out.t{1} == (3:1:8)') == 1);

% basic, non-edge case (iv input arguments)
one_ts_out = restrict(one_ts, rvec_1_iv);
assert(all(one_ts_out.t{1} == (3:1:8)') == 1);

% edge case -- note points equal to edges are INCLUDED
rvec_1 = [2, 8]; rvec_1_iv = iv(rvec_1); 

one_ts_out = restrict(one_ts, rvec_1(1), rvec_1(2));
assert(all(one_ts_out.t{1} == (2:1:8)') == 1);

one_ts_out = restrict(one_ts, rvec_1_iv);
assert(all(one_ts_out.t{1} == (2:1:8)') == 1);

%% preserve tsd data type
test_tsd = tsd();
test_tsd.tvec = 0:10;
test_tsd.data = zeros(size(test_tsd.tvec));

rvec = [0 10]; rvec_iv = iv(rvec);

test_tsd_out = restrict(test_tsd, rvec(1), rvec(2));
assert(CheckTSD(test_tsd_out) == 1);

test_tsd_out = restrict(test_tsd, rvec_iv);
assert(CheckTSD(test_tsd_out) == 1);

%% tsd functionality (1-D data)
test_tsd = tsd();
test_tsd.tvec = 0:10;
test_tsd.data = zeros(size(test_tsd.tvec));

rvec = [2.5 8.5]; rvec_iv = iv(rvec);
test_tsd_out = restrict(test_tsd, rvec(1), rvec(2));

assert(all(test_tsd_out.tvec == (3:1:8)) == 1); % tvec gets restricted correctly?
assert(sum(test_tsd_out.data) == 0); % data doesn't get corrupted?

%% tsd functionality (2-D data)
test_tsd = tsd();
test_tsd.tvec = 0:10;
test_tsd.data = zeros(5, length(test_tsd.tvec));

rvec = [2.5 8.5]; rvec_iv = iv(rvec);
test_tsd_out = restrict(test_tsd, rvec(1), rvec(2));

assert(CheckTSD(test_tsd_out) == 1);
assert(all(test_tsd_out.tvec == (3:1:8)) == 1); % tvec gets restricted correctly?
assert(sum(test_tsd_out.data(:)) == 0); % data doesn't get corrupted?
