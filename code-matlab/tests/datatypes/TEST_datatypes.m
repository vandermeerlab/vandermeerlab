function tests = TEST_datatypes
tests = functiontests(localfunctions);
end

function TestDataTypeInitialization(testCase)

test_iv = iv();
verifyTrue(testCase,CheckIV(test_iv));

test_ts = ts();
verifyTrue(testCase,CheckTS(test_ts));

test_tsd = tsd();
verifyTrue(testCase,CheckTSD(test_tsd));

end



function setupOnce(testCase)

end

function teardownOnce(testCase)

end
