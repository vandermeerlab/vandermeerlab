function BulkSpellCheck(cfg_in)
%BULKSPELLCHECK Make sure file names match rat session folder names /
%    verify correct folder contents.
%
%   Checks only for the presence of the session folder name (sessionID) in
%   the file or subdirectory name. Lists folder contents that do not match
%   the session folder name. Does not check spelling for any additional
%   prefixes or suffixes; does not recognize underscores in keys.m. Checks 
%   only 1 folder level deeper past the session folder. 
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
%
% Open and read the section "About this function" for more information
%
% ACarey May 2015, initial version 

%%

cfg_def.type = 'file';
cfg_def.directory = []; % can be a list 

cfg = ProcessConfig2(cfg_def,cfg_in);

%% helper functions

function findterminalfolder(cfg)
    dir_list = dir(pwd); % get me the folder info for all folders
    dir_list = dir_list(arrayfun(@(x) x.name(1), dir_list) ~= '.'); % removes ghost folders
 
    for iDir = 1:length(dir_list) 
        if dir_list(iDir).isdir 
            orig = pwd;
            cd([pwd,'\',dir_list(iDir).name])
            [~,sessionID,~] = fileparts(pwd);
            if  length(dir_list(iDir).name) == 15
                checkfiles(cfg,sessionID)
            else
                findterminalfolder(cfg);
            end
            cd(orig)
        end
    end
end

    function checkfiles(cfg,sessionID)
        %[~,sessionID,~] = fileparts(pwd);
        allListings = dir(pwd);
        allListings = allListings(arrayfun(@(x) x.name(1), allListings) ~= '.'); % removes ghost folders
        for iListing = 1:length(allListings)
            if ~allListings(iListing).isdir
                correct = strfind(allListings(iListing).name,sessionID);
                if isempty(correct)
                    disp([allListings(iListing).name, ' in ',sessionID])  
                end
            else
                backtoprev = pwd;
                cd(allListings(iListing).name)
                [~,folderName,~] = fileparts(pwd);
                allSubListings = dir(pwd);
                allSubListings = allSubListings(arrayfun(@(x) x.name(1), allSubListings) ~= '.'); % removes ghost folders
                for iSubListing = 1:length(allSubListings)
                    correct = strfind(allSubListings(iSubListing).name,sessionID);
                    if isempty(correct)
                        disp([allSubListings(iSubListing).name, ' in ',sessionID,'\',folderName])
                    end
                end
                cd(backtoprev)
            end
        end
        
    end

%% do the thing
originalFolder = pwd;
disp('BulkSpellCheck: Listing file names that do not match session folder names:')
if isempty(cfg.directory)
    findterminalfolder(cfg)
    
else % user has given a list of directories directly containing the files to be renamed
    for iDirectory = 1:length(cfg.directory)
        currentDirectory = cfg.directory{iDirectory};
        cd(currentDirectory)
        checkfiles(cfg)
    end

end
disp('BulkSpellCheck: finished checking.')
cd(originalFolder)
end

