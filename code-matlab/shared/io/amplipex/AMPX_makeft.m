function data_ft = AMPX_makeft(data)
% function data_ft = AMPX_makeft(data)
%
% creates fieldtrip compatible data structure from that returned by
% AMPX_loadData
%
% e.g. do
% channels_to_load = [10:13];
% data = AMPX_loadData(fname,channels_to_load,20);
% data_ft = AMPX_makeft(data);
%
% MvdM 2013-09-20

cfg = [];

data_ft = [];
data_ft.time = {data.tvec'};
data_ft.fsample = data.hdr.Fs;

for iC = length(data.channels):-1:1
   
    data_ft.trial{1}(iC,:) = data.channels{iC}';
    data_ft.label{iC} = num2str(data.labels(iC)); % may need to change this based on mapping
    
end

data_ft.sampleinfo = [1 length(data_ft.trial{1})];

data_ft.hdr.Fs = data.hdr.Fs;
data_ft.hdr.nchans = length(data.channels);
data_ft.hdr.nTrials = 1;
data_ft.hdr.nSamplesPre = 0;
data_ft.hdr.nSamples = length(data_ft.trial{1});
data_ft.hdr.chantype = {'unknown'};
data_ft.hdr.chanunit = {'unknown'};
data_ft.hdr.orig = data.hdr;

data_ft.cfg = cfg;
