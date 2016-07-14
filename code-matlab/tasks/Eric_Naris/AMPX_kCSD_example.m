function [csd_struct] = AMPX_kCSD_example(cfg_in, all_data_pre)
%% AMPX_kCSD_example: 
% takes a piece of data (idealy filtered and only an event of interest) and
% applies the kernal Current Source Density 

% get the input variables

%% Preamble 
cfg_def.save_path = '/Users/ecarmichael/Documents/Temp Data/CSD';
cfg_def.buffer = 0.5;

cfg = ProcessConfig2(cfg_in, cfg_def);

ExpKeys = all_data_pre.(cfg.session_name).lg.cycles.ExpKeys;
evts = all_data_pre.(cfg.session_name).lg.evts;

cfg.chan_to_plot = diag(ExpKeys.Probe_layout);
if strcmp(ExpKeys.ProbeType, 'A8x8')
    cfg.spacing = 0:.2:1.4;
end

%% extract the data for the example gamma event
fname = strrep(cfg.session_name, '_', '-');
cd(['D:\DATA\' fname(1:4) '\' fname(1:15) ])
fname = strrep(cfg.session_name, '-', '_');
[data, ~] = AMPX_Naris_preprocess([],fname,'pre');
data_remap_AMPX = AMPX_remap_sites(data, ExpKeys);


data_tsd = AMPX2tsd(data_remap_AMPX);
clear data
clear data_remap_AMPX

cfg_filter = [];
cfg_filter.f = [40 55];
cfg_filter.type = 'cheby1';
cfg_filter.order = 5;
event_data = FilterLFP(cfg_filter, data_tsd);

event_data_filt = restrict(event_data, evts.tstart(cfg.example)-cfg.buffer/2, evts.tend(cfg.example)+cfg.buffer/2);
event_data_not_filt = restrict(data_tsd, evts.tstart(cfg.example)-cfg.buffer/2, evts.tend(cfg.example)+cfg.buffer/2);


rm_chan = ExpKeys.BadChannels;
for iChan = rm_chan
    event_data_filt.data(iChan,:) = NaN;
end

%% plot the filter trace
c_ord = linspecer(8);
figure(1); 
subplot(211);
hold on
for iChan = 1:8
    plot(event_data_not_filt.tvec,event_data_not_filt.data(cfg.chan_to_plot(iChan),:), 'color', c_ord(iChan, :));
end
xlim([event_data_not_filt.tvec(1) event_data_not_filt.tvec(end)])
set(gca, 'ytick', [])

subplot(212);
hold on
for iChan = 1:8
    plot(event_data_filt.tvec,event_data_filt.data(cfg.chan_to_plot(iChan)), 'color', c_ord(iChan, :));
end
xlim([event_data_filt.tvec(1) event_data_filt.tvec(end)])
set(gca, 'ytick', [])

%% collect the event data and tvec
filt_data = event_data_filt.data;
filt_tvec = event_data_filt.tvec;
clear event_data_not_filt data_tsd event_data 
%% interpolate over missing channels
if exist('rm_chan','var') && ~isempty(rm_chan)
   for iS = 1:length(filt_data)
    temp_data = reshape(filt_data(:,iS), 8,8);
    temp_data = inpaintn(temp_data);
    filt_data(:,iS) = reshape(temp_data, 1,64);
   end
end


%% set up x and y coordinates for probe layout [NOTE: needs to match d.label when reshaped]
xspace = cfg.spacing; yspace = cfg.spacing;
xgrid = repmat(xspace',[1 8]); ygrid = repmat(yspace,[8 1]);

clear el_pos;
el_pos(:,1) = xgrid(:); % kCSD toolbox expects this format, see tutorial in doc folder
el_pos(:,2) = ygrid(:);
%%


nSamples = length(filt_tvec);
%for iS = nSamples:-1:1
loop = 1;
samples = 1:nSamples;
tic
for iS = 500:600% 1:length(samples);
    
    fprintf('Sample %d/%d...\n',iS,nSamples);
    
    this_sample = filt_data(:,iS);
    
    pots = this_sample(:); % kCSD expected format
    
    k = kcsd2d(el_pos, pots);
    
    X = k.X; Y = k.Y; CSD_est(iS,:,:) = k.CSD_est;
   
    CSD_1D(:,loop) = diag(squeeze(CSD_est(iS,:,:)));
    loop = loop+1;
end
t_sec = toc/60;
disp(['Elapsed time is ' num2str(t_sec) '  minutes'])
%% plot
rmpath('D:\Users\mvdmlab\My_Documents\GitHub\fieldtrip\external\stats')
subplot(2,2,1)
hold on
vline([filt_tvec(samples(1)), filt_tvec(samples(end))], {'r', 'r'}, {'Start', 'End'})
subplot(2,2,2)
imagesc(X(1,:),Y(:,1),squeeze(CSD_est(iS,:,:)));
colorbar
c_lim = get(gca, 'clim');
%% make a 2D plot using a specified layout. 
figure
subplot(2,1,1)
imagesc(CSD_1D);

set(gca, 'yTickLabel', {'DM', 'VL'}, 'yTick', [1 length(diag(squeeze(CSD_est(iS,:,:))))])
set(gca, 'xTickLabel', filt_tvec(500:600))
% colorbar

csd_struct.csd_1d = CSD_1D;
csd_struct.csd_est = CSD_est;
csd_struct.X = X;
csd_struct.Y = Y;
csd_struct.samples = samples;
csd_struct.cfg = cfg;
