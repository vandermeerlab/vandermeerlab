function pass_flag = tscmp(TS1,TS2)
%% tscmp Compare TS objects
%   TF = strcmp(TS1,TS2) compares the ts structs ts1 and ts2 and returns logical 1 (true)
%   if they are identical, and returns logical 0 (false) otherwise.
%
%   See also STRCMP
%
% youkitan 2016-03-22 initial version
% youkitan 2016-08-05 input check

%% check input
if ~CheckTS(TS1)
   error('ts1 is not a correctly formed ts.'); 
end

if ~CheckTS(TS2)
   error('ts2 is not a correctly formed ts.'); 
end

%% check overall structure
if length(TS1.t) ~= length(TS2.t);
    fprintf('There are unequal numbers of cells')
    pass_flag = false;
    return;
end

if length(TS1.label) ~= length(TS2.label);
    fprintf('There are unequal numbers of labels')
    pass_flag = false;
    return;
end

%% check content

if any(~cellfun(@isequal,TS1.t,TS2.t));
    fprintf('The spike times are not the same across all cells... \n')

    idx = find(~cellfun(@isequal,TS1.t,TS2.t));
    fprintf(['The following cells are not the same between files: ',num2str(idx),'\n'])
    
    pass_flag = false;
    return;
end
    
    
%% return true if all checks pass
pass_flag = true;

end