function fn = FindFile(globfn, varargin)

% fn = FindFile(globfn, parameters)
%
% Finds a single files that match a wildcard input globfn.
% Based on matlab's dir function.
% Can searche all directories under the current directory.
%
% INPUTS
%      globfn -- filename to search for (you can use '*',
%                but not '?', don't use directory names.
%
% OUTPUTS
%       fn -- string array of file found
%
% ERRORS
%       If zero or >1 file is found this produces an error.
%
% PARAMETERS
%       StartingDirectory -- default '.'
%       CheckSubdirs (1/0) -- default 0 -- note different from FindFiles
%
% ADR 1998
% version 1.0
%
% Status: PROMOTED (Release version) 

%-----------------
StartingDirectory = '.';
CheckSubdirs = 0;
extract_varargin;
fns = FindFiles(globfn, 'StartingDirectory', StartingDirectory, 'CheckSubdirs', CheckSubdirs);

if isempty(fns)
	error('FINDFILE:ZeroFiles', 'Zero files found.');
elseif length(fns) >1
	error('FINDFILE:TooManyFiles', 'More than one file found.');
else
    fn = fns{1};
end