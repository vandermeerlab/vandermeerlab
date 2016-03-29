function butter_check
% The ft butter function is different from the matlab one and will cause
% problems.  This will check to see if there is a butter function in one of
% the fieldtrip directories and remove it. 
%
%
% EC 08/01/16


path_butter = which('butter');
if isempty(findstr(path_butter, 'ield'))==0
rmpath(path_butter(1:end-8));
fprintf(['FieldTrip butter fcn removed from path \n'])
end