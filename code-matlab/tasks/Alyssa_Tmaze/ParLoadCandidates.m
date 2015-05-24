function evt = ParLoadCandidates
%ParLOADCANDIDATES Load candidates into BASE workspace 'transparently'
%for use in parfor loops
%   No need for any inputs; it checks the current directory automatically
%
% The candidates file must contain only one variable called 'evt' for this
% to work; when loaded, the variable appears in the workspace as 'evt'.
%
% Will not load anything if more than one file matching *candidates.mat is
% found in the current directory or subdirectories.
%
% A.Carey May 2015 

%%

%[~,name,~] = fileparts(pwd); 

fn = FindFiles('*candidates.mat');
if isempty(fn)
    evt = [];
    disp(['LoadCandidates: No files matching',' ''*candidates.mat'' ','were found in', [' ',pwd]])
    disp('candidates not loaded')
elseif length(fn) > 1
    evt = [];
    disp(['LoadCandidates: More than one file matching',' ''*candidates.mat'' ','was found in', [' ',pwd]])
    disp('candidates not loaded')
else
    evt = [];
    load(fn{1})
    assignin('base','evt',evt)
end

end

