function [data_out] = make_fake_data(cfg_in)
%% AMPX_fake_data: generates constant sine waves with a specified frequency, 
%  amplidue and phase offset.  Use for testings LFP analyes. Can output 
% tsd, AMPX, or FT formats using cfg.output
%
% Inputs
%   -cfg [struct]: contains all the parameters of the waves
%        cfg.nChan:number of channels (eg: 1:64 or [1 3 4 5]...)
%        cfg.Fs: sampling frequency
%        cfg.amp: amplitude
%        cfg.fs: frequency (constant at this point)
%        cfg.phase_off: phase offset from pervious channel
%        cfg.t_len: length of the wave/tvec (default: 10s)
%        cfg.output: fake data output format. Can be 'FT', 'TSD', or 'AMPX'

% Outputs
%   fake_data [struct]:
%         tvec
%         data
%         hdr
%         labels
%
%
% EC 02/10/2015
%
% Please feel free to expand since this is a function that is not called
% for any useful analysis...at least it never should be. 




%% make the wave for n channels with some offset
cfg_def.nChan = 1:64;
cfg_def.Fs = 2000; % in Hz
cfg_def.amp = 200;
cfg_def.amp_inc = 0;
cfg_def.freq = 52; % in Hz
cfg_def.phase_off = 0; % in degrees
cfg_def.t_len = 1; % time in seconds
cfg_def.output = 'AMPX';
cfg_def.plot = 'on';
cfg = ProcessConfig2(cfg_def, cfg_in);


c_ord = linspecer(length(cfg.nChan));

t = 0:1/cfg.Fs:cfg.t_len;

%% make the wave
wave = cell(length(cfg.nChan),1); 
labels = NaN*ones(1,length(cfg.nChan));

for iChan = 1:length(cfg.nChan)
    wave{iChan} = (cfg.amp+(cfg.amp_inc+iChan*5))*sin(2*pi*cfg.freq*t+(cfg.phase_off*iChan).*pi./180) ;
    labels(iChan) = cfg.nChan(iChan);
    if strcmp(cfg.plot, 'on')
    hold on
    plot(t, wave{iChan}, 'color', c_ord(iChan,:))
    end
end

%% basic data to be converted if needed (AMPX default)
        data.channels = wave;
        data.tvec =t;
        data.hdr.Fs = cfg.Fs;
        data.labels = labels;
%%
switch cfg.output
    
    case 'TSD'
        disp('Output TSD')
        data_out = AMPX2tsd(data);
%         fake_evts_Ids = 0:20:length(t);
%         fake_evts_s = fake_evts_Ids(2:end-1)-.05;
%         fake_evts_e = fake_evts_Ids(2:end-1)+.05;
%         fake_evts = [fake_evts_s; fake_evts_e]';
        
    case 'FT'
        disp('Output FT')
        data_out = AMPX_makeft(data);
        
    otherwise
        disp('Output AMPX')
        data_out = data;
        
end

