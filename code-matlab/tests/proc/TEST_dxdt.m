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

function TestType(testCase)

    % see if the type of the input variable remains the same after going
    % through the function
    tvec = 0:0.5:50;
    % create x array
    x = rand(1, length(tvec));
    start_type = class(x);
    dx = dxdt(tvec, x);
    end_type = class(dx);
    verifyEqual(testCase, start_type, end_type);

end

