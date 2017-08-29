function LoadCandidates(varargin)
%LOADCANDIDATES Load candidates into caller workspace
%   No need for any inputs; it checks the current directory automatically
%
% The candidates file must contain only one variable called 'evt' for this
% to work;when loaded, the variable appears in the workspace as 'evt'.
%
% Will not load anything if more than one file matching *candidates.mat is
% found in the current directory or subdirectories.
%
% A.Carey May 2015 
% MvdM 2017 modified to accept cfg.suffix input; implemented as varargin for
% compatibility

suffix = [];
switch nargin
    case 0
    case 1
        if isfield(varargin{1},'suffix')
            suffix = varargin{1}.suffix;
        else
            error('I couldn''t find any input fields I know how to handle.');
        end
    otherwise
        error('I don''t know what to do with more than one input argument.');
end

%[~,name,~] = fileparts(pwd); 

fstring = cat(2,'*-candidates',suffix,'.mat');
fn = FindFiles(fstring);
if isempty(fn)
    fprintf('LoadCandidates: no files matching %s were found in %s\n',fstring,pwd);
elseif length(fn) > 1
    fprintf('LoadCandidates: more than one file matching %s was found in %s\n',fstring,pwd);
    disp('candidates not loaded')
else
    load(fn{1})
    assignin('caller','evt',evt)
end

end

