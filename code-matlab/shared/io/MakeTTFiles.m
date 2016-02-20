function MakeTTFiles(cfg_in)
%MAKETTFILES Write .tt files to the current directory 
%   A .tt file contains all the spikes recorded by a single tetrode. Such 
%   a file may be useful for multiunit activity detection since the whole
%   tetrode file contains more information about spiking events than all
%   spike sorted units combined. Make sure to remove bad or noisy tetrodes 
%   from the file list ^_^
% 
%   .tt files can be loaded with LoadSpikes().
%
%   MAKETTFILES(cfg)
%
%   CONFIG OPTIONS
%     cfg.fc = {}; cell array containing full file names; this is the
%          output from FindFiles. If empty, MakeTTFiles will search for 
%          all .ntt files.
%     cfg.verbose = 1; 1 Talk to me, 0 don't
%
%   OUTPUT
%     files saved directly to current directory with identical names 
%     as the .ntt files, but with the extension as .tt
%
%   EXAMPLE USAGE:
%
%   CASE 1
%   % first, set your current directory to the session you want to work with
%   MakeTTFiles([]); % make tt file for all ntt in the current directory
%
%   CASE 2
%   % first, set your current directory to the session you want to work with
%   cfg = [];
%   cfg.fc = FindFiles('*.ntt');
%   cfg.fc(10) = []; % remove a bad tetrode from the file list
%   MakeTTFiles(cfg);
%
%   Sales pitch:
%     Ever look at your raster plot and wonder why you have so many fewer 
%     spikes than you saw when you were targeting cells? Well it's because 
%     you threw them all out when you were spike sorting! MakeTTFiles grabs
%     all the spikes recorded by a .ntt file and puts them in a format similar 
%     to a .t file. This way we have as much spiking data as possible when 
%     we are interested in detecting regions with lots of MUA. 
%     (More information == better?) 
% 
% aacarey Sept 2015
% parts modified from ADR's WriteTFiles (MClust function)

cfg_def.fc = {};
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% potential anger avoidance?: 
if ~isa(cfg.fc,'cell')
    error('MakeTTFiles: cfg.fc should be a cell array of filenames.');
end

% get files if none specified
if isempty(cfg.fc) 
    cfg.fc = FindFiles('*.ntt');
    if isempty(cfg.fc) 
        disp('MakeTTFiles: zero files found.') % fyi
        return
    end
end

if cfg.verbose; disp([mfun,': writing tt files...']); disp(' '); end

% do the thing
for iNTT = 1:length(cfg.fc) % for each ntt file
   % Load the timestamps (TS) [n x 1] double, but not the waveforms [n x 4 x 32] double where the '~' is (each spike "snapshot" taken by neuralynx consists of 32 samples)
   [TS,~] = LoadTT_NeuralynxNT(cfg.fc{iNTT}); 
   % convert timestamp units 
   %TS = TS./10^4; 
   
   % get the name of the .ntt file so it can be used for the .tt file
   [~,name,~] = fileparts(cfg.fc{iNTT});
   
   % make a new filename with .tt extension
   fn = [name,'.tt'];
   
   % begin to write the file named fn
   fp = fopen(fn, 'wb', 'b');
   
   % give it a header
   WriteTTHeader(fp);
   
   % finish writing the file
   fwrite(fp, TS, 'uint32');
   
   % close the file?
   fclose(fp);
end

disp(' '); if cfg.verbose; disp([mfun,': tt files written']); end

%~~~~~~~~~ helpers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% this nested function mimicks the header found in a .t file; modified from
% MClust's WriteHeader()
function WriteTTHeader(fp)
    fprintf(fp, '%%%%BEGINHEADER\n');
    fprintf(fp, '%% Program: matlab\n');
    fprintf(fp, '%% Date: %s\n', datestr(now));
    fprintf(fp, '%% Directory: %s\n', pwd);
    if ~isempty(getenv('COMPUTERNAME'))
        fprintf(fp, [ '%% Hostname: ', getenv('COMPUTERNAME'), '\n']);
    end
    if ~isempty(getenv('USERNAME'))
        fprintf(fp, [ '%% User: ', getenv('USERNAME'), '\n']);
    end
    fprintf(fp, '%% TT-file\n');
    fprintf(fp, '%% Timestamps for all spikes recorded by a tetrode\n');
    fprintf(fp, '%%%%ENDHEADER\n');
end
end

