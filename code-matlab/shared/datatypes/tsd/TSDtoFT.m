function data_ft = TSDtoFT(cfg_in,data)
% function data_ft = TSDtoFT(cfg_in,data_tsd)
%
% converts tsd into fieldtrip (ft) data structure
%
% cfg_def.mode = 'as-is'; % {'as-is','resample'}, defines how to deal with gaps in data
%
% MvdM 2014-11-12 initial version
% NOTE multiple channels not yet implemented!

cfg_def = [];
cfg_def.mode = 'as-is'; % {'as-is','resample'}, defines how to deal with gaps in data
cfg = ProcessConfig(cfg_def,cfg_in);

if ~CheckTSD(data)
   return; 
end

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

%
data_ft = [];
switch cfg.mode
    
    case 'as-is'
        nSamples = size(data.data,2);
        
        data_ft.trial{1}  = data.data;
        data_ft.time{1}   = data.tvec;
        
    case 'resample' % currently does interpolation of data, could be improved with options like inserting NaNs or zeros if no sample nearby
        
        data_ft.time{1} = data.tvec(1):1./Fs:data.tvec(end);
        data_ft.trial{1} = interp1(data.tvec,data.data,data_ft.time{1},'nearest');
                
        nSamples = length(data_ft.time{1});
    otherwise
        error('Unknown mode.');
end

data_ft.hdr.Fs = Fs;
data_ft.hdr.nSamples = nSamples;
data_ft.label   = data.label;
data_ft.sampleinfo = [1 data_ft.hdr.nSamples];