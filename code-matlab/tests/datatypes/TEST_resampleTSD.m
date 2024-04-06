function tests = TEST_datatypes
tests = functiontests(localfunctions);
end

function Test_resampleTSD_identical_rebinning(testCase)

cfg = []; cfg.method = 'rebinning';

% trivial case, don't mess up data if input times are equal to requested
% resampling times
new_tvec = (1:10)';
tsd_out = resampleTSD(cfg, testCase.TestData.tsd_in_1d, new_tvec);
verifyTrue(testCase, isequaln(testCase.TestData.tsd_in_1d.data, tsd_out.data));
end

function Test_resampleTSD_downsample_rebinning(testCase)
cfg = []; cfg.method = 'rebinning';

new_tvec = (1:2:10)'; % downsample
tsd_out = resampleTSD(cfg, testCase.TestData.tsd_in_1d, new_tvec);
verifyTrue(testCase, isequal(length(tsd_out.data), length(new_tvec)));

end

function Test_resampleTSD_identical_interp(testCase)

end


function setupOnce(testCase)

testCase.TestData.tsd_in_1d = tsd;
testCase.TestData.tsd_in_1d.tvec = (1:10)';
testCase.TestData.tsd_in_1d.data = rand(1,10);

testCase.TestData.tsd_in_2d = tsd;
testCase.TestData.tsd_in_2d.tvec = (1:10)';
testCase.TestData.tsd_in_2d.data = rand(2,10);

end

function teardownOnce(testCase)

end
