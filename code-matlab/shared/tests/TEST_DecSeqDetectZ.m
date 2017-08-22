function tests = DecSeqDetectZTest
% tests for DecSeqDetectZ() using MATLAB's unit testing framework

tests = functiontests(localfunctions);

end


function testFindSingleSequence(testCase)

% construct test sequence
nBins = 100;

testSeqData = 1:11; % bins that will be set to 1
testSeqT = 1:11;

P = tsd;
P.data = zeros(nBins,length(testSeqT));

for iT = 1:length(testSeqT)
    P.data(testSeqData(iT),testSeqT(iT)) = 1;
end
P.tvec = testSeqT;

% expected outcome
expected_iv = iv;
expected_iv.tstart = 1; expected_iv.tend = 11;

% run sequence detection
cfg = [];
seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Start times not equal');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'End times not equal');
end


function testFindSingleSequence(testCase)

% construct test sequence
nBins = 100;

testSeqData = 1:11; % bins that will be set to 1
testSeqT = 1:11;

P = tsd;
P.data = zeros(nBins,length(testSeqT));

for iT = 1:length(testSeqT)
    P.data(testSeqData(iT),testSeqT(iT)) = 1;
end
P.tvec = testSeqT;

% expected outcome
expected_iv = iv;
expected_iv.tstart = 1; expected_iv.tend = 11;

% run sequence detection
cfg = [];
seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Start times not equal');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'End times not equal');
end


function setupOnce(testCase)

fprintf('\nStarting tests...\n');
%warning('off','MATLAB:dispatcher:nameConflict'); % this doesn't get rid of
%path warning -- seems path is set before setupOnce is being called by
%runtests()

end

function teardownOnce(testCase)

fprintf('\nAll done!\n');

end