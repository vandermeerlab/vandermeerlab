function fd = getDataPath(cfg_in)
% GETDATAPATH Returns the base data path or a list of promoted folders.
%  fd = GETDATAPATH(cfg) Searches for specified promoted rat folders in the
%  base data directory. If no rats are specified, the base data directory is 
%  returned.
%  Note: to use this function, you need to edit it to include your own
%  computer's name and base data path.
%
%  INPUTS
%    cfg: config struct with fields controlling function behavior
%
%  OUTPUT
%    fd: cell array with found data folders
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
else
    machinename = getenv('HOSTNAME');
    filesep = '/';
end

% Choose base data path based on machine name.
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
    otherwise
        error('Unrecognized machine. Edit this function: add your machine name and base data path.')
end

fd = {};

curr_pwd = pwd;

if ~isempty(cfg.rats)
    for iRat = 1:length(cfg.rats)
        
        cd(base_fp);
        cd(cfg.rats{iRat});
        
        temp_fd = dir;
        temp_fd = temp_fd(3:end); % remove . and ..
        temp_fd = temp_fd([temp_fd.isdir]); % only consider dirs
        
        for iFD = 1:length(temp_fd)
            
            cd(temp_fd(iFD).name);
                                  
            % accept
            fd = cat(1,fd,pwd);
            
            cd .. % return to rat folder
            
        end % of session folders
    end % of rats
else
    fd = base_fp;
end

cd(curr_pwd) % return to starting folder

