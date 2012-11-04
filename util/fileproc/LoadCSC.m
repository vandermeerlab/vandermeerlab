function [csc,csc_info] = LoadCSC(fname,varargin)
% function [csc,csc_info] = LoadCSC(fname,varargin)
%
% MvdM 12

TimeConvFactor = 10^-6; % from nlx units to seconds
VoltageConvFactor = 10^6; % from volts to microvolts
extract_varargin;

% load data
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

% extract information from header
csc_info = readCSCHeader(Header);

% unwrap samples
csc_data = reshape(Samples,[size(Samples,1)*size(Samples,2) 1]);

% apply conversion
csc_data = csc_data.*VoltageConvFactor.*csc_info.ADBitVolts;

% create matching timestamps; idea is to make a matrix to match Samples
% original and then add dt's within each 512 sample packet
csc_timestamps = repmat(Timestamps,[size(Samples,1) 1]).*TimeConvFactor;

dtvec = (0:size(Samples,1)-1)*(1/csc_info.SamplingFrequency);
dtmat = repmat(dtvec',[1 size(Samples,2)]);

csc_timestamps = csc_timestamps+dtmat;
csc_timestamps = reshape(csc_timestamps,[size(csc_timestamps,1)*size(csc_timestamps,2) 1]);

% done
csc = tsd(csc_timestamps,csc_data);



function csc_info = readCSCHeader(Header)

csc_info = [];
for hline = 1:length(Header)
   
    line = strtrim(Header{hline});
    
    if isempty(line) | ~strcmp(line(1),'-') % not an informative line, skip
        continue;
    end
    
    a = regexp(line(2:end),'(?<key>\w+)\s+(?<val>\S+)','names');
    
    % deal with characters not allowed by MATLAB struct
    if strcmp(a.key,'DspFilterDelay_µs')
        a.key = 'DspFilterDelay_us';
    end
    
    csc_info = setfield(csc_info,a.key,a.val);
    
    % convert to double if possible
    if ~isnan(str2double(a.val))
        csc_info = setfield(csc_info,a.key,str2double(a.val));
    end
    
end