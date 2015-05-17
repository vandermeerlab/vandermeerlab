function GatherFiles(cfg_in)
%GATHERFILES Copy or move multiple files with a similar name
%     All files matching a wildcard input in the provided search directories
%     or current subdirectories (default) can be moved into a single directory. 
%     You will be prompted to respond [y/n] using keyboard input before the
%     files are moved or copied.
%     Originally intended for gathering analysis image files from different rat
%     sessions into a single folder for viewing. 
%   GatherFiles(cfg_in)
%       cfg_def.directory = ''; % list of directories; if empty, checks all 
%                subfolders in current directory
%       cfg_def.globfn = ''; % wildcard input, the globfn of FindFiles()
%       cfg_def.destination = ''; % the directory that the files will be
%                moved/copied to
%       cfg_def.mode = 'copy'; 'move' 
%       cfg_def.prompt = 1; % if 0, the action specified in cfg.mode
%                will be performed without prompting you for permission
%
% ACarey May 2015

%% use

%cd('D:\data\promoted')

%cfg.globfn = '*EventPower.png';
%cfg.destination = 'D:\data\temp';
%cfg.mode = 'move';

%GatherFiles(cfg);

%%
cfg_def.directory = ''; % list of directories; if empty, checks all subfolders in current directory
cfg_def.globfn = '';
cfg_def.destination = '';
cfg_def.mode = 'copy';
cfg_def.prompt = 1;

cfg = ProcessConfig2(cfg_def,cfg_in);

if isempty(cfg.destination)
    error('cfg.destination must be specified')
end

originalFolder = pwd;

if isempty(cfg.directory)
    fns = FindFiles(cfg.globfn);
else
    for iDirectory = 1:length(cfg.directory)
        currentDirectory = cfg.directory{iDirectory};
        cd(currentDirectory)
        fns = FindFiles(cfg.globfn);
    end
end

if isempty(fns)
    disp('No files found')
else

    disp(' ')
    disp(['GatherFiles found the following',[' ',num2str(length(fns))],' file(s) to',' ',cfg.mode,' to ',cfg.destination,':'])
    disp(' ')
    arrayfun(@(x) disp(x), fns)
        if cfg.prompt
            %disp('Any keypress other than Y will not move the file(s)')
            prompt = ['Would you like to',[' ',cfg.mode],' the file(s)? Y/N (Enter):'];
            str = input(prompt,'s');
        else
            str = 'y';
        end
    if strcmp(str,'y') || strcmp(str,'Y')
        switch cfg.mode
            case 'move'
                action = 'moved';
                for iFile = 1:length(fns)
                    movefile(fns{iFile},cfg.destination)
                end
            case 'copy'
                action = 'copied';
                for iFile = 1:length(fns)
                    copyfile(fns{iFile},cfg.destination);
                end
        end
        
        disp(['Files',' ',action,' to',' ',cfg.destination])
    else
        disp('Files not moved or copied')
    end
    cd(originalFolder)
end

end