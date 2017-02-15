function ls_out = lsdep(mfun)
%% LSDEP List Dependencies 
%   LSDEP(mfun) lists the main dependencies for the function(s) named in mfun and returns
%   a cell array with all function paths.
%      
%   lsdep(S), checks for all function calls in the function named by string S. 
%    
%   lsdep(C), checks for all unique function calls in the functions named in cell array C.
%
% youkitan - 2016-12-09

if ischar(mfun)
    [f_all,~] = matlab.codetools.requiredFilesAndProducts(mfun);
    
elseif iscell(mfun)
    f_all = {};
    for iC = 1:length(mfun)
        [this_f,~] = matlab.codetools.requiredFilesAndProducts(mfun{iC});
        f_all = [f_all this_f];
    end
    
else
    error('incorrect input')
end

for i = 1:length(f_all)
    disp(f_all{i})
end
    
ls_out = f_all;
end