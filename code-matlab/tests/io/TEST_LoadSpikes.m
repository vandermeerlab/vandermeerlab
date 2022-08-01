function tests = TEST_LoadSpikes
tests = functiontests(localfunctions);
end

function BasicDataLoadingTest(testCase)

testCase.TestData.S1 = LoadSpikes([]);

verifyTrue(testCase,CheckTS(testCase.TestData.S1));
verifyEqual(testCase,length(testCase.TestData.S1.t),92);

cfg = []; cfg.load_questionable_cells = 1;
testCase.TestData.S2 = LoadSpikes(cfg);

verifyTrue(testCase,CheckTS(testCase.TestData.S2));
verifyEqual(testCase,length(testCase.TestData.S2.t),128);
end



function setupOnce(testCase)
testCase.TestData.origPath = pwd;

%fprintf('TEST_LoadSpikes is called from %s.\n',pwd);
%ws = getenv('WORKSPACE');
%fprintf('getenv workspace is %s.\n',ws);

testFolder1 = fullfile(pwd,'testdata','R050-2014-04-02');
cd(testFolder1)
end

function teardownOnce(testCase)
cd(testCase.TestData.origPath);
end
