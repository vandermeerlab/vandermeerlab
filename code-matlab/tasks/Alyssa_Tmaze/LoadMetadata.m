function LoadMetadata
%LOADMETADATA Load metadata struct into caller workspace
%   No need for any inputs; it checks the current directory automatically
%
% The metadata file must contain only one variable called 'metadata' for this
% to work as expected. The variable metadata may contain any number of
% additional fields. 
%
% The folder and metadata file must both be named according to the standard
% RXXX-XXXX-XX-XX (ratID-year-month-day). 
%
% This function will fail if, for example, your metadata file is called
% R050-2014-04-03-metadata.mat but the folder that holds it is
% called R050-2014-04-03_recording. Make sure your folders and data are
% named according to the standard (as outlined above). 
%
% A.Carey Feb, 2015; written for T-maze-specific metadata file

%%

[~,name,~] = fileparts(pwd); 

fn = FindFiles([name,'-metadata.mat']);
load(fn{1})
assignin('caller','metadata',metadata)

end

