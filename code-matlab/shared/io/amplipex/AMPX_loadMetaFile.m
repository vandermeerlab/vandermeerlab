function hdr = AMPX_loadMetaFile(fname,varargin)
% function hdr = AMPX_loadMetaFiler(fname,varargin)
%
% load amplipex .meta file
%
% INPUTS:
% fname: full filename of .meta file
%
% MvdM 2013-03-27

extract_varargin;

keys = {'Amplitude range max', ...
        'Amplitude range min', ...
        'File length (Amplipex clock (sec)', ...
        'File size (bytes)', ...
        'TimeStamp of the end of recording (computer clock - ms)', ...
        'TimeStamp of the start of recording (computer clock - ms)', ...
        'Number of recorded channels', ...
        'Recording start date', ...
        'Recording start time', ...
        'Sampling rate', ...
        'Sha1 code for file'};

values = {'range1_volts', ...
          'range2_volts', ...
          'filelength_sec', ...
          'filesize_bytes', ...
          'end_timestamp', ...
          'start_timestamp', ...
          'nChannels', ...
          'date', ...
          'time', ...
          'Fs', ...
          'SHA1'};

fid = fopen(fname);

if fid == -1
   error('Could not open file %s',fname); 
end

hdr = [];

iL = 1;
while 1
    
   line = fgetl(fid);
   if line == -1
       break;
   end
   
   hdr.raw{iL} = line;
   iL = iL + 1;
   
   tokens = regexp(line,'(.*)=(.*)','tokens');
   tokens = strtrim(tokens{1});
   
   if ~isnan(str2double(tokens{2})) % check if convertible to numeric
       tokens{2} = str2double(tokens{2});
   end
   
   for iK = 1:length(keys)
       
       idx = strmatch(tokens{1},keys); % check if this key is on the list
       
       if ~isempty(idx) % on list, replace
           hdr = setfield(hdr,values{idx},tokens{2});
       else % not on list, use literally
           hdr = setfield(hdr,tokens{1},tokens{2});
       end
      
   end
   
end

fclose(fid);