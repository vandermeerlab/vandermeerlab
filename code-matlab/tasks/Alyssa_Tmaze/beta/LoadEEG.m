function LFP = LoadEEG(cfg_in,channel)
%LOADEEG Loads Buzsaki lab .eeg file containing local field potential data
%   LFP = LOADEEG(cfg_in,channel)
%
%    INPUTS
%       cfg: config struct with fields controlling function behavior
%       channel: channel number you want to load
%
%    OUTPUTS
%       LFP: tsd struct
%
%    CONFIG OPTIONS
%       cfg.VoltageConvFactor = 10^-7;
%       cfg.verbose = 1;
%
%  Make sure to addpath for Buzsaki lab's code.
%
% aacarey Dec 2015 from MvdM code

%cfg_def.channel = [];
%cfg_def.TimeConvFactor = 1;
cfg_def.VoltageConvFactor = 10^-7;
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose; fprintf('%s: Loading channel data...\n',mfun); end

xml = LoadXml(FindFile('*.xml'));

[Data, OrigIndex] = LoadBinary(FindFile('*.eeg'),channel);

tvec = OrigIndex .* (1/xml.lfpSampleRate);
Data = Data.*cfg.VoltageConvFactor;

LFP = tsd(tvec,Data);

MakeHeader

LFP = History(LFP,mfun,cfg);


    function MakeHeader
       hdr.FileType = 'EEG';
       hdr.SamplingFrequency = xml.lfpSampleRate;
       % more here?
       
       LFP.cfg.hdr = {hdr};
    end
end

