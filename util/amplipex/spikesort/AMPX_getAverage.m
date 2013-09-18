function av = AMPX_getAverage(fname,varargin)
% function av = AMPX_getAverage(fname,iChan,decimate_factor)
%
% returns average signal computed from specified channels
%
% NOTE spike filtering is off by default
%
% MvdM 13

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

spk_filter = 0;

% status
fprintf('AMPX_getAverage(): loading %d channels...\n',length(iChan));
tic;

% amount of bytes to skip after reading each sample
szINT16 = 2;  % sizeof(int16)=2
skipBytes = (nbChan-1)*szINT16;
nSamples = Inf; % number of samples to read, Inf for all

fid = fopen(fname,'rb');

% loop through channels to load
for iC = length(iChan):-1:1
        
    offset = (iChan(iC)-1)*szINT16; % compute how many bytes to skip
    fseek(fid, offset, 'bof'); % move reading head to first sample of desired signal
    
    current_channel = fread(fid,nSamples,'*int16',skipBytes);
    current_channel = double(current_channel);
    
    if decimate_factor > 1
       
        current_channel = decimate(current_channel,decimate_factor);
        
    end
    
    if spk_filter
        current_channel = filter_for_spikes(current_channel);
    end
    
    if iC == length(iChan) % first channel loaded
        av = current_channel;
    else
        av = av + current_channel;
    end
    
end

av = av ./ length(iChan);

fclose(fid);
fprintf('AMPX_getAverage(): that took %.2f seconds.\n',toc);