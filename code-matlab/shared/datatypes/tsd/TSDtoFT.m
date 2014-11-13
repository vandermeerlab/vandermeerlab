function data_ft = TSDtoFT(cfg_in,data)
% function data_ft = TSDtoFT(cfg_in,data_tsd)
%
% converts tsd into fieldtrip (ft) data structure
%
% MvdM 2014-11-12 initial version

cfg_def = [];
cfg = ProcessConfig2(cfg_def,cfg_in);

if ~CheckTSD(data)
   return; 
end

nSamples = size(data.data,2);

dts = unique(diff(data.tvec));
if length(dts) > 1
   fprintf('\nTSDtoFT.m: WARNING: tvec diffs are not constant, cannot determine Fs.');
   % could approximate with a median if matches average closely enough
   Fs = 1./median(dts);
   fprintf('\nTSDtoFT.m: Fs %.2f estimated.\n',Fs);
else
   Fs = 1./dts;
   fprintf('\nTSDtoFT.m: Fs %.2f detected.\n',Fs);
end

data_ft           = [];
data_ft.trial{1}  = data.data;
data_ft.time{1}   = data.tvec;

data_ft.hdr.Fs = Fs;
data_ft.hdr.nSamples = nSamples;

data_ft.label   = data.label;
data_ft.sampleinfo = [1 data_ft.hdr.nSamples];

