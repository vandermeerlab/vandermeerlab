function fd = getTmazeDataPath(cfg_in)
% function fd = getTmazeDataPath(cfg_in)
%
% get list of folders with data to analyze
%
% assumes that data is organized as AllDataFolder > RatFolder >
% SessionFolder, i.e. C:\data\R050\R050-2014-04-03 etc..
%
% CONFIGS
%
% cfg_def.rats = {'R042','R044','R050'};
% cfg_def.requireMetadata = 1;
%
% OUTPUT
%
% fd: cell array with found data folders
%
% MvdM 2015

cfg_def.rats = {'R042','R044','R050'};
cfg_def.requireMetadata = 1;

cfg = ProcessConfig2(cfg_def,cfg_in);

if ispc
    machinename = getenv('COMPUTERNAME');
    filesep = '\';
else
    machinename = getenv('HOSTNAME');
    filesep = '/';
end

switch machinename
    
    case 'ISIDRO'
        base_fp = 'C:\data\';
    case 'EQUINOX'
        base_fp = 'D:\vandermeerlab\data\';
    case 'ATHENA'
        base_fp = 'D:\vandermeerlab\data\';
end

fd = {};

curr_pwd = pwd;

for iRat = 1:length(cfg.rats)
   
    cd(base_fp);
    cd(cfg.rats{iRat});
    
    temp_fd = dir; 
    temp_fd = temp_fd(3:end); % remove . and ..
    temp_fd = temp_fd([temp_fd.isdir]); % only consider dirs
    
    for iFD = 1:length(temp_fd)
        
        cd(temp_fd(iFD).name);
        
        if cfg.requireMetadata
           m = FindFiles('*metadata.mat');
           if isempty(m)
               cd ..
               continue;
           end
        end
        
        % accept
        fd = cat(1,fd,pwd);
        
        cd .. % return to rat folder
        
    end % of session folders
    
end

cd(curr_pwd) % return to starting folder

