function pass_flag = tccmp(TC1,TC2)
%% TCCMP Compare Tuning Curve objects
%   TF = tccmp(TS1,TS2) compares the tc structs tc1 and tc2 and returns logical 1 (true)
%   if they are identical, and returns logical 0 (false) otherwise.
%
%   See also STRCMP
%
% youkitan 2016-12-02 initial version

%% check input
mfun = mfilename;

if ~CheckTC(TC1,mfun)
   error('tc1 is not a correctly formed tc.'); 
end

if ~CheckTC(TC2,mfun)
   error('tc2 is not a correctly formed tc.'); 
end

%% check overall structure
if ~isequal(size(TC1.tc),size(TC2.tc))
    fprintf('The sizes of the tuning curves are unequal.\n')
    pass_flag = 0;
    return
end

if ~isequal(size(TC1.pos_idx),size(TC2.pos_idx))
    fprintf('The sizes of the binned tuning variables are unequal.\n')
    pass_flag = 0;
    return
end

if isfield(TC1,'usr') ~= isfield(TC2,'usr')
    fprintf('<<< WARNING: Only one of the tc objects has usr fields! >>>\n')    
end

% IMPLEMENT MORE TESTS...

%% check content
if ~isequal(TC1.tc,TC2.tc)
    fprintf('The contens of the tuning curves are not equal.\n')
    pass_flag = 0;
    return
end

if ~isequal(TC1.spk_hist,TC2.spk_hist)
    fprintf('The contents of the spike counts are not equal.\n')
    pass_flag = 0;
    return
end

if ~isequal(TC1.occ_hist,TC2.occ_hist)
    fprintf('The contents of the occupancies are not equal.\n')
    pass_flag = 0;
    return
end

if ~isequal(TC1.pos_idx,TC2.pos_idx)
    fprintf('The contents of the binned tuning variables are not equal.\n')
    pass_flag = 0;
    return
end
    
% IMPLEMENT MORE TESTS...

%% return true if all checks pass
pass_flag = true;

end