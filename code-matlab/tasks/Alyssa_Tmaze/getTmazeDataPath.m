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
% cfg_def.rats = {'R042','R044','R050','R064'};
% cfg_def.requireMetadata = 1;
% cfg_def.requireCandidates = 0;
% cfg_def.requireEvents = 0;
% cfg_def.verbose = 1;
% cfg_def.userpath = ''; if specified, uses this path instead of default
%
% OUTPUT
%
% fd: cell array with found data folders
%
% MvdM 2015
% youkitan 2016-11-22 edit: added user input for path

cfg_def.rats = {'R042','R044','R050','R064'};
cfg_def.requireMetadata = 1;
cfg_def.requireCandidates = 0;
cfg_def.requireEvents = 0;
cfg_def.verbose = 1;
cfg_def.userpath = '';

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if ispc
    machinename = getenv('COMPUTERNAME');
    filesep = '\';
elseif ismac
    machinename = getenv('USER');
    filesep = '/';
else
    machinename = getenv('HOSTNAME');
    filesep = '/';
end

% overide default if user specifies path for data folder
if ~isempty(cfg.userpath)
    machinename = 'USERDEFINED';
end

switch machinename
    
    case 'ISIDRO'
        base_fp = 'C:\data\';
    case {'EQUINOX','BERGKAMP'}
        base_fp = 'D:\data\';
    case 'MVDMLAB-ATHENA'
        base_fp = 'D:\vandermeerlab\';
    case {'MVDMLAB-EUROPA','DIONYSUS'}
        base_fp = 'D:\data\promoted\';
    case 'CALLISTO'
        base_fp = 'E:\data\promoted\';
    case 'USERDEFINED'
        base_fp = cfg.userpath;
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
        
        if cfg.requireCandidates
           m = FindFiles('*candidates.mat');
           if isempty(m)
               cd ..
               continue;
           end
        end
        
        if cfg.requireEvents
            m = FindFiles('*.nev');
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

