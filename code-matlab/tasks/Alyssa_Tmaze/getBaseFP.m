function base_fp = getBaseFP
% function base_fp = getBaseFP
%
% get the location of source data folders
%
%
% OUTPUT
%
% base_fp: the directory of source data folders on your computer. Note: if
% your computer is not recognized by getBaseFP you need to modify the
% function and tell it where you store your promoted data files.
%
% aacarey Nov 2015, from getTmazeDataPath (MvdM)
% aacarey 02 MAR 2019, edit to give additional instructions to new users

if ispc
    machinename = getenv('COMPUTERNAME');
    %filesep = '\';
else
    machinename = getenv('HOSTNAME');
    %filesep = '/';
end

switch machinename
    % You can get your PC's machine name by typing "getenv('COMPUTERNAME')"
    % in the command window. Then, add your computer's name to the cases
    % below, or create a new case for your data path. The data path is
    % whichever folder you've saved vandermeerlab rat data (R042, R044,
    % R050 and R064 folders).
    case 'ISIDRO'
        base_fp = 'C:\data\';
    case {'EQUINOX','BERGKAMP','VYSERITHUS'}
        base_fp = 'D:\data\';
    case 'MVDMLAB-ATHENA'
        base_fp = 'D:\vandermeerlab\';
    case {'MVDMLAB-EUROPA','DIONYSUS'}
        base_fp = 'D:\data\promoted\';
    case 'CALLISTO'
        base_fp = 'E:\data\promoted\';
    otherwise
        error('Unknown computer. Edit this function to add your data path and run the analysis.')
end

