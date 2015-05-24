function metadata = ParLoadMetadata
%PARLOADMETADATA Load metadata struct into BASE workspace 'transparently'
%for use in parfor loops
%   No need for any inputs; it checks the current directory automatically
%
% The metadata file must contain only one variable called 'metadata' for this
% to work as expected. The variable metadata may contain any number of
% additional fields.
%
% ACarey May 2015 

%%

%[~,name,~] = fileparts(pwd); 

fn = FindFiles('*metadata.mat');
if isempty(fn)
    metadata = [];
    disp(['ParLoadMetadata: No files matching',' ''*metadata.mat'' ','were found in', [' ',pwd]])
    disp('metadata not loaded')
elseif length(fn) > 1
    metadata = [];
    disp(['ParLoadMetadata: More than one file matching',' ''*metadata.mat'' ','was found in', [' ',pwd]])
    disp('metadata not loaded')
else
    metadata = [];
    load(fn{1})
    assignin('base','metadata',metadata)
end

end

