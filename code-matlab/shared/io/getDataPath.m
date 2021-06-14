function [fd, fd_extra] = getDataPath(cfg_in)
% GETDATAPATH Returns the base data path or a list of promoted folders.
%  fd = GETDATAPATH(cfg) Searches for specified promoted rat folders in the
%  base data directory. If no rats are specified, the base data directory is 
%  returned.
%
%  Note: to use this function, you need to edit it to include your own
%  computer's name and base data path.
%
%  INPUTS
%    cfg: config struct with fields controlling function behavior
%
%  OUTPUTS
%    fd: cell array with found data folders
%    fd_extra: struct with additional information:
%     .rat_ID: rat IDs (cell array)
%     .rat_ID_num: numerical rat IDs (double)
%     .fd_date: dates (datetime array) ***NOTE, ASSUMES THAT DATA FOLDERS
%       END WITH THEIR DATE***
%     .fd_date_num: numerical dates, starting with the first date as 1
%
%  CONFIG OPTIONS
%    cfg_def.rats = {}; Cell array of strings specifying which rats to get
%                       data paths for. Ex: {'R042','R050'}. Returns a list 
%                       of all session directories for these rats. If cfg.rats 
%                       is empty {}, returns the base data path.
%
%          * Assumes that data is organized as AllDataFolder > RatFolder >
%            SessionFolder, i.e. C:\data\R050\R050-2014-04-03 etc..
%
% MvdM 2015
% aacarey edit Jan 2016

cfg_def.rats = {};
cfg_def.verbose = 1;

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

% Choose base data path based on machine name.
switch machinename
    
    case {'ISIDRO','MVDMLAB-PERSEUS','ODYSSEUS'}
        base_fp = 'C:\data\';
    case {'EQUINOX','BERGKAMP'}
        base_fp = 'D:\data\adrlab';
    case 'MVDMLAB-ATHENA'
        base_fp = 'D:\vandermeerlab\';
    case {'MVDMLAB-EUROPA','DIONYSUS'}
        base_fp = 'D:\data\promoted\';
    case 'CALLISTO'
        base_fp = 'E:\data\promoted\';
    case {'MVDMLAB-PHOBOS'}
        base_fp = 'D:\ADRLabData';
    case {'manishm'}
        base_fp = '/Users/manishm/Work/vanDerMeerLab/ADRLabData';
    otherwise
        error('Unrecognized machine. Edit this function: add your machine name and base data path.')
end

fd = {};
ratID = {}; ratID_num = [];
fd_date = []; fd_date_num = [];

curr_pwd = pwd;

if ~isempty(cfg.rats)
    for iRat = 1:length(cfg.rats)
        
        % collect data for this rat
        cd(base_fp);
        cd(cfg.rats{iRat});
        
        temp_fd = dir;
        temp_fd = temp_fd(3:end); % remove . and ..
        temp_fd = temp_fd([temp_fd.isdir]); % only consider dirs
        
        this_fd = {};
        this_ratID = {}; this_ratID_num = [];
        this_fd_date = [];
        
        for iFD = 1:length(temp_fd)
            
            cd(temp_fd(iFD).name);
                                  
            % accept
            this_pwd = pwd;
            this_fd = cat(1,this_fd,this_pwd);
            this_ratID = cat(1,this_ratID,cfg.rats{iRat});
            this_ratID_num = cat(1,this_ratID_num,iRat);
            
            % date
            this_fd_date = cat(1,this_fd_date,datetime(this_pwd(end-9:end))); % note, assumes date is at end of filename!
            
            cd .. % return to rat folder
            
        end % of session folders
        
        % sort data for this rat, and add to master lists
        [this_fd_date, sort_idx] = sort(this_fd_date);
        this_fd = this_fd(sort_idx);
        this_fd_date_num = days(this_fd_date - this_fd_date(1)) + 1; % days counting from first session date as 1
        
        fd = cat(1,fd,this_fd);
        ratID = cat(1,ratID,this_ratID);
        ratID_num = cat(1,ratID_num,this_ratID_num);
        fd_date = cat(1,fd_date,this_fd_date);
        fd_date_num = cat(1,fd_date_num,this_fd_date_num);
        
    end % of rats
else
    fd = base_fp;
end

fd_extra.rat_ID = ratID;
fd_extra.ratID_num = ratID_num;
fd_extra.fd_date = fd_date;
fd_extra.fd_date_num = fd_date_num;

cd(curr_pwd) % return to starting folder