function LoadMetadata
%LOADMETADATA Load metadata struct into caller workspace
%   No need for any inputs; it checks the current directory automatically
%
% The metadata file must contain only one variable called 'metadata' for this
% to work as expected. The variable metadata may contain any number of
% additional fields.  
%
% aacarey Feb 2015
%  --edit May 2015 

%%

fn = FindFiles('*metadata.mat');
if isempty(fn)
    disp(['LoadMetadata: No files matching ''*metadata.mat'' were found in ',pwd])
elseif length(fn) > 1
    disp(['LoadMetadata: More than one file matching ''*metadata.mat'' was found in ',pwd])
    disp('metadata not loaded')
else
    load(fn{1})
    assignin('caller','metadata',metadata)
end

end

