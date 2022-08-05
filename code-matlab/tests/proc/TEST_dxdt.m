function tests = TEST_dxdt()
    tests = functiontests(localfunctions);
    
end
% function passflag = TEST_dxdt()
%
% unit tests for dxdt()

function TestOutputSize(testCase)

    % create time vector
    tvec = 0:0.5:50;
    % create x array
    x = rand(1, length(tvec));

    dx = dxdt(tvec, x);
    verifyEqual(testCase,length(dx),length(x));

end

