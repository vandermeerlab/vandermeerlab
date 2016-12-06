function pass_flag = tccmp(TC1,TC2)
%% TCCMP Compare Tuning Curve objects
%   TF = tccmp(TS1,TS2) compares the tc structs tc1 and tc2 and returns logical 1 (true)
%   if they are identical, and returns logical 0 (false) otherwise.
%
%   See also STRCMP
%
% youkitan 2016-12-02 initial version

%% check input
if ~CheckTC(TC1)
   error('tc1 is not a correctly formed tc.'); 
end

if ~CheckTC(TC2)
   error('tc2 is not a correctly formed tc.'); 
end

%% check overall structure


%% check content

    
    
%% return true if all checks pass
pass_flag = true;

end