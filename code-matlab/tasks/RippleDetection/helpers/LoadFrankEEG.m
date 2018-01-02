function LFP = LoadFrankEEG(cfg,filename)
%LOADFRANKEEG Loads EEG data provided by Frank and colleages (from CRCNS)
%   LFP = LoadFrankEEG(cfg,filename)
%
%   CFG OPTIONS
%     cfg.verbose = 1 (default)
%          1 - Display text in command window
%          0 - Do not display text in command window
%
%   INPUTS
%         cfg - Config struct with fields controlling function behaviour
%    filename - String with extentsion specifying the the EEG data you want
%               to load (e.g. 'boneeg07-1-13.mat'); only loads one file at 
%               a time.
%
%   OUTPUTS
%         LFP - local field potential (EEG) data presented in a TSD
%               datatype format
%
% aacarey Dec 2017 initial version
% see also: LoadCSC

cfg_def.verbose = 1;
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);

% Make sure user has specified a file that exists in the directory
if isempty(FindFiles(filename))
    error('File does not exist in directory. Check spelling.');
end

% Get base file name
[~,name,~] = fileparts(filename);

% Don't load the wrong type of data
if ~strcmp(name(4:6),'eeg')
    error('The file is not EEG data')
end

% Talk to me
if cfg.verbose; fprintf('%s: loading %s...\n',mfun,filename); end

% Load the data
eeg = load(filename,'eeg');

% Get indices into data array
recordingDay = str2double(name(7:8));
epoch = str2double(name(10));
tetrode = str2double(name(12:13));

% Get voltage data
eegStruct = eeg.eeg{1,recordingDay}{1,epoch}{1,tetrode};
data = eegStruct.data;

% Construct tvec
stepSize = 1/eegStruct.samprate;
tvec = eegStruct.starttime:stepSize:((length(data)-1)/eegStruct.samprate + eegStruct.starttime);

% Make output TSD with label (name)
LFP = tsd(tvec',data',name);

% Add header info and cfg history
LFP.cfg.hdr{1,1}.SamplingFrequency = eegStruct.samprate;
LFP = History(LFP,mfun,cfg);

% Make sure TSD is well formed
if ~CheckTSD(LFP)
    error('Failed to load EEG data as a proper TSD datatype')
end

end

