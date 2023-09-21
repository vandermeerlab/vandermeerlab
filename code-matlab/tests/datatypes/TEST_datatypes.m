function tests = TEST_datatypes
tests = functiontests(localfunctions);
end

function TestDataTypeInitialization(testCase)

test_iv = iv();
verifyTrue(testCase, CheckIV(test_iv));

test_ts = ts();
verifyTrue(testCase, CheckTS(test_ts));

test_tsd = tsd();
verifyTrue(testCase, CheckTSD(test_tsd));

end

function TestTSFunctions(testCase)

test_ts = ts();
test_ts.label{1} = 'a'; test_ts.label{2} = 'b';
test_ts.t{1} = [1 2]; test_ts.t{2} = [3 4];

cfg = []; cfg.mode = 'next';
[t, idx] = FindTSTime(cfg, test_ts, 2.5);
verifyTrue(testCase, t == 3); verifyTrue(testCase, idx == 2);

cfg.mode = 'prev';
[t, idx] = FindTSTime(cfg, test_ts, 2.5);
verifyTrue(testCase, t == 2); verifyTrue(testCase, idx == 1);

end


function setupOnce(testCase)

end

function teardownOnce(testCase)

end
