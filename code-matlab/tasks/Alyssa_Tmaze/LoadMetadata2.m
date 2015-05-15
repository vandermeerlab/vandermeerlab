function tf = LoadMetadata2
%LOADMETADATA2 Load metadata struct into caller workspace
%   No need for any inputs; it checks the current directory automatically
%
%   tf = 0 if file not found, and 1 if file found.
%        if tf = 1 then metadata has been loaded into caller workspace
%
% The metadata file must contain only one variable called 'metadata' for this
% to work as expected
%
%  Used in scripts that generate metadata fields (need to check if metadata 
%  exists already)
%
% A.Carey Feb 2015
% --edit May 2015

%% example usage

%loaded = LoadMetadata2;

% if ~loaded % then it doesn't exist yet, so make a metadata struct
%     metadata.SWRtimes = SWRtimes;
%     metadata.SWRfreqs = SWRfreqs;
% else % it does, so add a new field
%     metadata.SWRtimes = SWRtimes;
%     metadata.SWRfreqs = SWRfreqs;
% end

% [~,name,~] = fileparts(pwd); % pwd is your current folder, we just want its namepart

% savename = strcat(name,'-metadata.mat'); 
% save(savename,'metadata'); 

%%

%[~,name,~] = fileparts(pwd); 

fn = FindFiles('*metadata.mat');

if isempty(fn)
    tf = 0;
elseif length(fn) > 1
    disp(['LoadMetadata2: More than one file matching',' ''*metadata.mat'' ','was found in', [' ',pwd]])
    disp('metadata not loaded')
else
    tf = 1;
    load(fn{1})
    assignin('caller','metadata',metadata)
end

