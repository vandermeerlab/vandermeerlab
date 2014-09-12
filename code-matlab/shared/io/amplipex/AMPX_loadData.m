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
% NOTE AMPX .dat file format is a_1 b_1 c_1 ... a_n b_n c_n
% where a, b, c are signal channels and 1...n are samples (int16)
%
% MvdM 2013

% set some constants
szINT16 = 2;  % sizeof(int16)=2
rgINT16 = (2^16)./2;

if ~exist(fname,'file')
    error('File not found')
end

% first load .meta file
meta_fname = strrep(fname,'dat','meta');
data.hdr = AMPX_loadMetaFile(meta_fname);

nbChan = data.hdr.nChannels;

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
data.hdr.Fs = data.hdr.Fs ./ decimate_factor;

% status
fprintf('AMPX_loadData(): loading %d channels...\n',length(iChan));
tic;

% amount of bytes to skip after reading each sample
skipBytes = (nbChan-1)*szINT16;
nSamples = Inf; % number of samples to read, Inf for all

fid = fopen(fname,'rb');

% loop through channels to load
for iC = length(iChan):-1:1
        
    offset = (iChan(iC)-1)*szINT16; % compute how many bytes to skip
    fseek(fid, offset, 'bof'); % move reading head to first sample of desired signal
    
    data.channels{iC} = fread(fid,nSamples,'*int16',skipBytes); % this is now a signed 16bit integer
    
    data.labels(iC) = iChan(iC);
    
    % convert to microvolts
    data.channels{iC} = double(data.channels{iC})./rgINT16; % convert to fraction of full range
    data.channels{iC} = data.channels{iC}.*data.hdr.range1_volts*10^6; % convert to microvolts
    data.channels{iC} = data.channels{iC}./data.hdr.Gain; % correct for amplifier gain
    
    % check if we need to decimate
    if decimate_factor > 1
       
        data.channels{iC} = decimate(data.channels{iC},decimate_factor);
        
    end
    
    % remove DC component (note this could also be done in the raw file..)
    data.channels{iC} = data.channels{iC}-mean(data.channels{iC});
       
end

% construct tvec
data.tvec = 0:1./data.hdr.Fs:data.hdr.filelength_sec;
data.tvec = data.tvec(1:end-1)';

fclose(fid);

fprintf('AMPX_loadData(): that took %.2f seconds.\n',toc);
