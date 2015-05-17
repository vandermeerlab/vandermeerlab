function BulkRename(cfg_in)
%BULKRENAME For those 2 times in a row you messed up 18 filenames and didn't
% want to have to rename them by hand a 7th time
%
%   By default, assumes that your session folders have been named according to the
%   15-character standard RXXX-XXXX-XX-XX (ratID-year-month-day), and that
%   your filenames share these 15 characters at least. You can specify
%   additional/different nameparts in the config options 
%   *** the full, correct name of the file must be described in the config
%   options (when added to the RXXX-XXXX-XX-XX part)
%
%   BulkRename(cfg_in)
%      cfg_def.type = 'file'; 
%      cfg_def.directory = []; % directories to look into, can be a list 
%            of fd's directly containing the files to be renamed; if empty,
%            searches all subdirectories in current directory 
%      cfg_def.subfolder = ''; % look in a subfolder in session folder
%            ex: 'files', to look in D:\data\promoted\R042\R042-2013-08-21\files
%      cfg_def.oldprefix = ''; % any characters appearing before the standard
%      cfg_def.oldsuffix = ''; % any characters appearing after the
%                              standard, can include the extension as well
%      cfg_def.oldID = ''; % if it's different from the standard
%      cfg_def.oldext = ''; % the old extension, if relevant
%      cfg_def.newID = ''; % if you would like to change the 15 standard
%                              characters
%      cfg_def.newprefix = ''; 
%      cfg_def.newsuffix = ''; 
%      cfg.newext = ''; % change the extension
%      cfg_def.verbose = 1; % tell me what you're renaming; 0, don't tell me
%
% Open and read the section "About this function" for more information
%
% ACarey May 2015, initial version 

%% About this function 

%%% When and how to use this function:

% Suppose you have a script that generates a file for each session for each
% rat, and you accidentally named the file incorrectly for nSessions, where 
% n happens to be a burdensomely large number:
%   R042-2013-08-16-metadata_text.mat_oops   ...etc
%   R044-2013-12-16-metadata_text.mat_oops   ...etc
%   R050-2014-03-29-metadata_text.mat_oops   ...etc
% 
% You can use this function to automatically rename these files, such that
% the suffix is replaced:
% cfg.oldsuffix ='_text.mat_oops';
% cfg.newsuffix = '_test';
% cfg.newext = '.mat';

% Or you can do it with fewer config fields:
% cfg.oldsuffix ='_text.mat_oops';
% cfg.newsuffix = '_test.mat';

%%% How this function works

% If the directory to look in is unspecified, function looks in all
% subfolders until it finds what it consider a terminal folder: a folder 
% with a 15-character-long name, which it *assumes* has the format
% RXXX-XXXX-XX-XX. Then it concatenates all the user's descriptions of the
% file names, looks for a file with that name, and renames it to what the
% user has specified. 

%%

cfg_def.type = 'file';
cfg_def.directory = []; % can be a list 
cfg_def.subfolder = '';
cfg_def.oldprefix = ''; % any characters appearing before the standard
cfg_def.oldsuffix = ''; % any characters appearing after the standard
cfg_def.oldext = ''; 
cfg_def.oldID = '';
cfg_def.newID = '';
cfg_def.newprefix = ''; 
cfg_def.newsuffix = '';
cfg_def.newext = '';
cfg_def.verbose = 1; % talk to me

cfg = ProcessConfig2(cfg_def,cfg_in);

%% helper functions

function findterminalfolder(cfg)
    dir_list = dir(pwd); % get me the folder info for all folders
    dir_list = dir_list(arrayfun(@(x) x.name(1), dir_list) ~= '.'); % removes ghost folders
 
    for iDir = 1:length(dir_list) 
        if dir_list(iDir).isdir 
            orig = pwd;
            cd([pwd,'\',dir_list(iDir).name])
            if  length(dir_list(iDir).name) == 15
                [~,sessionID,~] = fileparts(pwd);
                if ~isempty(cfg.subfolder)
                    subfolder = ['\',cfg.subfolder];
                    cd([pwd,subfolder])
                end
                renamefile(cfg,sessionID)
            else
                findterminalfolder(cfg);
            end
            cd(orig)
        end
    end
end

    function renamefile(cfg,sessionID)
        %[~,sessionID,~] = fileparts(pwd);
        if ~isempty(cfg.newID)
            newfilename = [cfg.newprefix,cfg.newID,cfg.newsuffix,cfg.newext];
        else
            newfilename = [cfg.newprefix,sessionID,cfg.newsuffix,cfg.newext];
        end
        
        filename = [cfg.oldprefix,sessionID,cfg.oldsuffix,cfg.oldext];
        fn = FindFiles(filename);
        if isempty(fn)
            disp(strcat('No files called',[' ',filename,],' were found in ',[' ',pwd]))
        elseif length(fn) > 1
            % maybe this is impossible but w/e:
            error(strcat('More than one file called',[' ',filename,],' was found in',[' ',pwd]))
        else             
            
            
            if cfg.verbose
                if ~isempty(cfg.oldID)
                    oldfilename = [cfg.oldprefix,cfg.oldID,cfg.oldsuffix,cfg.oldext];
                else
                    oldfilename = [cfg.oldprefix,sessionID,cfg.oldsuffix,cfg.oldext];
                end
                
                disp(strcat('Renaming ',[' ',oldfilename],' to ',[' ',newfilename],' in ',[' ',pwd]))
            end
            movefile(filename,newfilename)

        end
        
    end

%% do the thing
originalFolder = pwd;

if isempty(cfg.directory)
    findterminalfolder(cfg)
    
else % user has given a list of directories directly containing the files to be renamed
    for iDirectory = 1:length(cfg.directory)
        currentDirectory = cfg.directory{iDirectory};
        cd(currentDirectory)
        renamefile(cfg)
    end

end

cd(originalFolder)
end

