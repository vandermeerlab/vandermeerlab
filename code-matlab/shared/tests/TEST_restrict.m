function passflag = TEST_restrict(cfg_in)
% function passflag = TEST_restrict(cfg_in)
%
% unit tests for restrict()
%
% MvdM 2015-11-30 basic initial version


cfg_def = [];
cfg_def.verbose = 1;

cfg = ProcessConfig(cfg_def,cfg_in);


passflag = 1; % set to 1 by default, make 0 if any test fails

% preliminaries: report which version of restrict() we are using

fprintf('\n*** UNIT TESTS FOR RESTRICT ***\n');


fp = which('restrict');
fprintf('\nTesting %s\n\n',fp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BLOCK 1: TESTS FOR PRESERVING INPUT FORMAT CORRECTLY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct some testing data - TS
empty_ts = ts;

one_ts_data = (0:1:10)';
one_ts = ts; one_ts.t{1} = one_ts_data; one_ts.label{1} = 'label1';

rvec = [0 10]; rvec_iv = iv(rvec);

% test - TS
empty_ts_out = restrict(empty_ts,rvec(1),rvec(2));
empty_ts_out_flag = CheckTS(empty_ts_out);
switch empty_ts_out_flag
    case 1
        fprintf('PASSED: empty TS, raw restrict argument\n');
    case 0
        fprintf('FAILED: empty TS, raw restrict argument\n');
        passflag = 0;
end

empty_ts_out = restrict(empty_ts,rvec_iv);
empty_ts_out_flag = CheckTS(empty_ts_out);
switch empty_ts_out_flag
    case 1
        fprintf('PASSED: empty TS, iv restrict argument\n');
    case 0
        fprintf('FAILED: empty TS, iv restrict argument\n');
        passflag = 0;
end

% construct some testing data - TSD

% construct some testing data - IV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BLOCK 2: TESTS FOR CORRECT PROCESSING OF EDGE CASES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BLOCK 2: SPEED TESTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%