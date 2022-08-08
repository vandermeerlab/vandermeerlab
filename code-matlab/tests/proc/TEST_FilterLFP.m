function tests = TEST_FilterLFP
tests = functiontests(localfunctions);
end

function BasicFilterLFPTest(testCase)
% make fake lfp
nPoints = 10000;
fTSD = tsd(1:nPoints, rand(1, nPoints));

fTSD.cfg.hdr{1}.SamplingFrequency = 2000;

out = FilterLFP([], fTSD);

verifyTrue(testCase, CheckTSD(out)); % verify that CheckTSD gives us a proper TSD
verifyEqual(testCase, length(fTSD.data), length(out.data));  % verify that the lenght of the output is equal to the length of the input

end