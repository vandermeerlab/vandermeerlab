function LoadMetadata
%LOADMETADATA Load metadata struct into caller workspace
%   No need for any inputs; it checks the current directory automatically
%
% The metadata file must contain only one variable called 'metadata' for this
% to work as expected
%
% A.Carey Feb, 2015; written for T-maze-specific metadata file

%%

[~,name,~] = fileparts(pwd); 

fn = FindFiles([name,'-metadata.mat']);
load(fn{1})
assignin('caller','metadata',metadata)

end

