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

if ispc
    machinename = getenv('COMPUTERNAME');
    %filesep = '\';
else
    machinename = getenv('HOSTNAME');
    %filesep = '/';
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
end

