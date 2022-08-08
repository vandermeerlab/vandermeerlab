function tests = DecSeqDetectZTest
% tests for DecSeqDetectZ() using MATLAB's unit testing framework

tests = functiontests(localfunctions);

end


function testFindDefinedSequences(testCase)

% construct single test sequence; easiest case
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

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Single sequence start time incorrect');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'Single sequence end time incorrect');

% single NaN case - no skip
P.data(:,6) = NaN;
expected_iv.tstart = [1; 7]; expected_iv.tend = [5; 11];

cfg = []; seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Single noskip NaN sequence start times incorrect');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'Single noskip NaN sequence end times incorrect');

% skip
expected_iv.tstart = 1; expected_iv.tend = 11;

cfg = []; cfg.nMaxNanSkipSequential = 1; seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Single skip NaN sequence start time incorrect');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'Single skip NaN sequence end time incorrect');

% single NaN case followed by too big jump
P.data(:,7) = 0; P.data(50,7) = 1; 

expected_iv.tstart = [1; 8]; expected_iv.tend = [5; 11];

cfg = []; cfg.nMaxNanSkipSequential = 1;  seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Single skipjump NaN sequence start times incorrect');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'Single skipjump NaN sequence end times incorrect');

% double NaN case
P.data(:,7) = NaN;
expected_iv.tstart = [1; 8]; expected_iv.tend = [5; 11];

cfg = []; cfg.nMaxNanSkipSequential = 1;  seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Double skip NaN sequence start times incorrect');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'Double skip NaN sequence end times incorrect');

% check that double skip works
expected_iv.tstart = 1; expected_iv.tend = 11;

cfg = []; cfg.nMaxNanSkipSequential = 2;  seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Double 2skip NaN sequence start times incorrect');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'Double 2skip NaN sequence end times incorrect');

% double skip followed by jump
P.data(:,8) = 0; P.data(50,8) = 1; 
expected_iv.tstart = [1; 9]; expected_iv.tend = [5; 11];

cfg = []; cfg.nMaxNanSkipSequential = 2;  seq_iv = DecSeqDetectZ(cfg,P);

verifyEqual(testCase,seq_iv.tstart,expected_iv.tstart,'Double 2skipjump NaN sequence start times incorrect');
verifyEqual(testCase,seq_iv.tend,expected_iv.tend,'Double 2skipjump NaN sequence end times incorrect');

end


function testFindRandomSequences(testCase)

% construct test sequence
nBins = 100;
nPoints = 1e5;

rng('default');
testSeqData = randi(nBins,[1 nPoints]); % for each column, find bin that will be set to 1
testSeqIdx = ((0:nPoints-1)*nBins); % create index vector to set zero matrix elements to 1
testSeqIdx = testSeqIdx + testSeqData;

testSeqT = 1:nPoints;

P = tsd;
P.data = zeros(nBins,nPoints);
P.data(testSeqIdx) = 1;
P.tvec = testSeqT;

% run sequence detection -- first without NaN skips selected
cfg = [];
cfg.nMaxNanSkipSequential = 0;
cfg.minLength = 3;
seq_iv_nogap = DecSeqDetectZ(cfg,P); nSeq_nogap = length(seq_iv_nogap.tstart);
verifyGreaterThan(testCase,nSeq_nogap,0,'Number of detected sequences not greater than zero!');

% now with NaN skips selected
cfg.nMaxNanSkipSequential = 1;
seq_iv_gap1 = DecSeqDetectZ(cfg,P); nSeq_gap1 = length(seq_iv_gap1.tstart);

% #seqs for above 2 runs should be the same, because no NaNs in the data!
verifyEqual(testCase,nSeq_nogap,nSeq_gap1,'Gap and no-gap seq counts on non-NaN data not equal');

% now introduce some NaNs
nanIdx = randi(nPoints,[1 floor(nPoints./10)]); % find column idxs that will be set to NaN
P.data(:,nanIdx) = NaN;

seq_iv_gap1nan = DecSeqDetectZ(cfg,P); nSeq_gap1nan = length(seq_iv_gap1nan.tstart);

% check if any starts or ends are NaN
seq_start = P.data(1,seq_iv_gap1nan.tstart);
verifyFalse(testCase,any(isnan(seq_start)),'At least one 1-skip sequence starts with NaN');
seq_end = P.data(1,seq_iv_gap1nan.tend);
verifyFalse(testCase,any(isnan(seq_end)),'At least one 1-skip sequence ends with NaN');

% run no gaps allowed on NaN data
cfg.nMaxNanSkipSequential = 0;
seq_iv_gap0nan = DecSeqDetectZ(cfg,P); nSeq_gap0nan = length(seq_iv_gap0nan.tstart);
verifyLessThan(testCase,nSeq_gap0nan,nSeq_gap1nan,'Seq counts on NaN data not less for no-gap-allowed detection');

% check if any starts or ends are NaN
seq_start = P.data(1,seq_iv_gap0nan.tend);
verifyFalse(testCase,any(isnan(seq_start)),'At least one 0-skip sequence starts with NaN');
seq_end = P.data(1,seq_iv_gap0nan.tend);
verifyFalse(testCase,any(isnan(seq_end)),'At least one 0-skip sequence ends with NaN');

% verify lengths
len = seq_iv_nogap.tend-seq_iv_nogap.tstart + 1;
verifyEqual(testCase,cfg.minLength,min(len),sprintf('Smallest 0-skip seq (%d) detected on non-NaN random data not equal to cfg.minLength',min(len)));
len = seq_iv_gap1.tend-seq_iv_gap1.tstart + 1;
verifyEqual(testCase,cfg.minLength,min(len),sprintf('Smallest 1-skip seq (%d) detected non-NaN random data not equal to cfg.minLength',min(len)));
len = seq_iv_gap1nan.tend-seq_iv_gap1nan.tstart + 1;
verifyEqual(testCase,cfg.minLength,min(len),sprintf('Smallest 1-skip seq (%d) detected on 10%% NaN random data not equal to cfg.minLength',min(len)));
len = seq_iv_gap0nan.tend-seq_iv_gap0nan.tstart + 1;
verifyEqual(testCase,cfg.minLength,min(len),sprintf('Smallest 0-skip seq (%d) detected on 10%% NaN random data not equal to cfg.minLength',min(len)));

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