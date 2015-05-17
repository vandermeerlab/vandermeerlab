function BulkDelete(cfg_in)
%BULKDELETE Because normal people make mistakes sometimes
%  
%   By default, assumes that your session folders have been named according to the
%   15-character standard RXXX-XXXX-XX-XX (ratID-year-month-day), and that
%   your filenames share these 15 characters at least. You can specify
%   additional/different nameparts in the config options 
%   *** the full, correct name of the file must be described in the config
%   options (when added to the RXXX-XXXX-XX-XX part)
%
%   BulkDelete(cfg_in)
%      cfg_def.type = 'file'; 
%      cfg_def.directory = []; % directories to look into, can be a list 
%            of fd's directly containing the files to be renamed; if empty,
%            searches all subdirectories in current directory 
%      cfg_def.subfolder = ''; % look in a subfolder in session folder
%            ex: 'files', to look in D:\data\promoted\R042\R042-2013-08-21\files
%      cfg_def.prefix = ''; % any characters appearing before the standard
%      cfg_def.suffix = ''; % any characters appearing after the
%                              standard, can include the extension as well
%      cfg_def.ext = ''; % extension, if relevant
%
% Open and read the section "About this function" for more information
%
% ACarey May 2015, initial version 

%% About this function 

%%% When and how to use this function:

% Use this function to delete unneeded files from your session folders:
%   R042-2013-08-16-metadata_test.mat   ...etc
%   R044-2013-12-16-metadata_test.mat   ...etc
%   R050-2014-03-29-metadata_test.mat   ...etc
% 
% You can specify these files using config options (it automatically finds
% the sessionID from the session folders)
% cfg.suffix ='_test';
% cfg.ext = '.mat';

% Or you can do it with fewer config fields:
% cfg.suffix ='_test.mat';

%%% How this function works

% If the directory to look in is unspecified, function looks in all
% subfolders until it finds what it considers a terminal folder: a folder 
% with a 15-character-long name, which it *assumes* has the format
% RXXX-XXXX-XX-XX. Then it concatenates all the user's descriptions of the
% file names, looks for a file with that name, and deletes it. 

%%

cfg_def.type = 'file';
cfg_def.directory = []; % can be a list 
cfg_def.subfolder = '';
cfg_def.prefix = ''; % any characters appearing before the standard
cfg_def.suffix = ''; % any characters appearing after the standard
cfg_def.ext = ''; 
cfg_def.verbose = 1; % talk to me

cfg = ProcessConfig2(cfg_def,cfg_in);

%% helper functions

    function deletelist = FindFilesToDelete(cfg,deletelist)
        dir_list = dir(pwd); % get me the folder info for all folders
        dir_list = dir_list(arrayfun(@(x) x.name(1), dir_list) ~= '.'); % removes 'ghost' folder things
        
        for iDir = 1:length(dir_list)
            if dir_list(iDir).isdir
                orig = pwd;
                cd([pwd,'\',dir_list(iDir).name])
                if  length(dir_list(iDir).name) == 15
                    [~,sessionID,~] = fileparts(pwd);
                    filename = [cfg.prefix,sessionID,cfg.suffix,cfg.ext];
                    disp(cfg.subfolder)
                    if ~isempty(cfg.subfolder)
                        subfolder = ['\',cfg.subfolder];
                        cd([pwd,subfolder])
                    end
                    %fn = FindFiles(filename);
                    deletelist = [deletelist FindFiles(filename)];
                else
                    deletelist = FindFilesToDelete(cfg,deletelist);
                end
                cd(orig)
            end
        end
    end

%% do the thing
originalFolder = pwd;

if isempty(cfg.directory)
    deletelist = [];
    deletelist = FindFilesToDelete(cfg,deletelist);
    
else % user has given a list of directories directly containing the files to be renamed
    for iDirectory = 1:length(cfg.directory)
        currentDirectory = cfg.directory{iDirectory};
        cd(currentDirectory)
        fname = [cfg.prefix,sessionid,cfg.suffix,cfg.ext];
        %fn = FindFiles(filename);
        deletelist = [deletelist FindFiles(fname)];
    end
end

if isempty(deletelist)
    disp('0 files found')
else
    disp(' ')
    disp(['BulkDelete found the following',[' ',num2str(length(deletelist))],' file(s) to delete:'])
    %disp(deletelist)
    arrayfun(@(x) disp(x), deletelist)
    %disp('Any keypress other than Y will not delete the file(s)')
    prompt = 'Would you like to delete the file(s)? Y/N (Enter):';
    str = input(prompt,'s');
    if strcmp(str,'y') || strcmp(str,'Y')
        for iFile = 1:length(deletelist)
                delete(deletelist{iFile})
        end
        disp('Files deleted')
    else
        disp('Files not deleted ^_^')
    end
cd(originalFolder)
end

end

