function csc_tsd = LoadCSC_Axona(cfg_in)
% function csc_tsd = LoadCSC_Axona(cfg)
%
% loads Neuralynx .ncs files, handling missing samples correctly
%
% INPUTS: only cfg
%
% cfg.fc: cell array containing filenames to load
%   if no file_list field is specified, loads all *.Ncs files in current dir
% cfg.verbose = 1; Allow or suppress displaying of command window text
%
% OUTPUTS:
%
% csc_tsd: CSC data tsd struct
%
% MvdM 2018-07-26
% based on Sturla Molden's getegf.m

cfg_def.fc = {};
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % this takes fields from cfg_in and puts them into cfg

if isempty(cfg.fc) % no filelist provided, load everything
    
    cfg.fc = FindFiles('*.egf');
    
else
    
    if ~isa(cfg.fc,'cell')
        error('LoadCSC_Axona: cfg.fc should be a cell array.');
    end
    
end

if isempty(FindFiles(cfg.fc{1}))
    error('File does not exist in directory. Check spelling.');
end

fc = sort(cfg.fc);
nFiles = length(fc);

% track sample counts for different files, throw error if not equal
sample_count_tvec = nan(nFiles,1);
sample_count_data = nan(nFiles,1);

csc_tsd = tsd;

if cfg.verbose; fprintf('%s: Loading %d file(s)...\n',mfun,nFiles); end

for iF = 1:nFiles
    
    fname = fc{iF};
    
    % load raw data
    [data, Fs] = getegf(fname);

    % construct tvec
    tvec = (0:length(data)-1)'/Fs;

    % check if the data is the same length for each channel.  
    if iF > 1 && length(data) ~= length(csc_tsd.data(iF-1,:))
        error('Data lengths differ across channels.');
    end

    % update cfg info
    csc_tsd.cfg.hdr{iF}.Fs = Fs;
    
    % done, add to tsd
    csc_tsd.tvec = tvec;
    csc_tsd.data(iF,:) = data;
    
    [~,fn,fe] = fileparts(fname);
    csc_tsd.label{iF} = cat(2,fn,fe);

end

% add sessionID
[~,csc_tsd.cfg.SessionID,~] = fileparts(pwd);

% housekeeping
csc_tsd.cfg.Fs = Fs;
csc_tsd.cfg.history.mfun = cat(1,csc_tsd.cfg.history.mfun,mfun);
csc_tsd.cfg.history.cfg = cat(1,csc_tsd.cfg.history.cfg,{cfg});

function [eeg, Fs] = getegf(filename)
% Read EGF file
%
% [EEG,Fs] = geteef(datafile);
%
% You can create a time vector for plotting your EEG by taking 
% t = (0:length(EEG)-1)'/Fs
%
% This routine requires the Signal Processing Toolbox! 
%
%
% Sturla Molden <sturla@molden.net>
% Centre for the Biology of Memory
% Norwegian University of Science and Technology
% http://www.cbm.ntnu.no
% 
% Copyright (C) 2003  Centre for the Biology of Memory, NTNU
% All Rights Reserved
%
% This M-file is released under Q Public License v 1.0,
% a copy of which should accompany this file.

fid = fopen(filename,'r','ieee-le');

if (fid == -1)
    error(sprintf('Cannot open file %s',filename));
end

for ii = 1:8
    string = fgetl(fid);
end
Fs = sscanf(string,'%*s %u %*s');
for ii = 1:2
    string = fgetl(fid);
end
nsamp = sscanf(string,'%*s %u');
fseek(fid,10,0);
eeg = fread(fid,nsamp,'int16');
fclose(fid);
if (Fs == 960)
  	eeg = decimate(eeg,2,100,'FIR');
    Fs = 480;
elseif (Fs == 4800)
    eeg = decimate(eeg,8,100,'FIR');
    Fs = 600;
else
   fclose(fid);
   error('Unknown sampling rate');
end