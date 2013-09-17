function data = AMPX_loadData(fname,varargin)
% function data = AMPX_loadData(fname,iChan,decimate_factor)
%
% load raw Amplipex data
% NOTE: consider running AMPX_RemoveDCfromDat.m first!
%
% INPUTS:
%
% fname: .dat file to read
% iChan: array of channels to read (if omitted, read everything)
% decimate_factor: decimate data by specified factor (if omitted, do not decimate)
% NOTE iChan needs to be specified for decimate_factor to work
%
% OUTPUTS:
%
% data: [1 x 1 struct] with fields as follows
% .channels: {nChannels x 1} cell array with data [nSamples x 1 int16] in each cell
% .labels: [nChannels x 1] channel IDs matching .channels
% .hdr: [1 x 1 struct] header read from .meta file
% .tvec: [nSamples x 1 double] time axis
%
% MvdM 2013

if ~exist(fname,'file')
    error('File not found')
end

% first load .meta file
meta_fname = strrep(fname,'dat','meta');
data.hdr = AMPX_loadMetaFile(meta_fname);

nbChan = str2double(data.hdr.nChannels); % this is clumsy, should be fixed in AMPX_loadMetaFile

% figure out what to load
if nargin == 1 % read all channels
    iChan = 1:nbChan;
else
    iChan = varargin{1};
end

if nargin == 3
    decimate_factor = varargin{2};
else
    decimate_factor = 1;
end
data.hdr.Fs = str2double(data.hdr.Fs) ./ decimate_factor;

m = memmapfile(fname,'Format','int16','writable',false);

% loop through channels to load
for iC = length(iChan):-1:1
        
    data.channels{iC} = m.Data(iChan(iC):nbChan:end);
    data.labels(iC) = iChan(iC);
    
    % check if we need to decimate
    if decimate_factor > 1
       
        data.channels{iC} = decimate(double(data.channels{iC}),decimate_factor);
        
    end
       
end

% construct tvec
data.tvec = 0:1./data.hdr.Fs:str2double(data.hdr.filelength_sec);
