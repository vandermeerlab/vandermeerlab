function [csd_struct] = AMPX_kCSD(cfg_in, filt_data)
%% AMPX_kCSD: 
% takes a piece of data (idealy filtered and only an event of interest) and
% applies the kernal Current Source Density 

%% Preamble 
LoadExpKeys
cfg_def.save_path = '/Users/ecarmichael/Documents/Temp Data/CSD';
cfg_def.chan_to_plot = [38 43 29 22 10 2 56 57];
cfg_def.buffer = 0;
if strcmp(ExpKeys.ProbeType, 'A8x8')
cfg_def.spacing = 0:.2:1.4; 

end

cfg = ProcessConfig2(cfg_in, cfg_def);

rm_chan = ExpKeys.BadChannels;
for ichan = rm_chan
    filt_data.data(ichan,:) = NaN;
end

figure(1); subplot(221);
plot(filt_data.tvec,filt_data.data(1,:));
%% interpolate over missing channels
if exist('rm_chan','var') && ~isempty(rm_chan)
   for iS = 1:length(filt_data.data)
    temp_data = reshape(filt_data.data(:,iS), 8,8);
    temp_data = inpaintn(temp_data);
    filt_data.data(:,iS) = reshape(temp_data, 1,64);
   end
end


%% set up x and y coordinates for probe layout [NOTE: needs to match d.label when reshaped]
xspace = cfg.spacing; yspace = cfg.spacing;
xgrid = repmat(xspace',[1 8]); ygrid = repmat(yspace,[8 1]);

clear el_pos;
el_pos(:,1) = xgrid(:); % kCSD toolbox expects this format, see tutorial in doc folder
el_pos(:,2) = ygrid(:);
%%


nSamples = length(filt_data.tvec);
%for iS = nSamples:-1:1
loop = 1;
samples = 1:nSamples;
tic
for iS = 1:length(samples);
    
    fprintf('Sample %d/%d...\n',iS,nSamples);
    
    this_sample = filt_data.data(:,iS);
    
    pots = this_sample(:); % kCSD expected format
    
    k = kcsd2d(el_pos, pots);
    
    X = k.X; Y = k.Y; CSD_est(iS,:,:) = k.CSD_est;
   
    CSD_1D(:,loop) = diag(squeeze(CSD_est(iS,:,:)));
    loop = loop+1;
end
toc
%% plot
subplot(2,2,1)
hold on
vline([filt_data.tvec(samples(1)), filt_data.tvec(samples(end))], {'r', 'r'}, {'Start', 'End'})
subplot(2,2,2)
imagesc(X(1,:),Y(:,1),squeeze(CSD_est(iS,:,:)));
colorbar
c_lim = get(gca, 'clim');
%% make a 2D plot using a specified layout. 

subplot(2,2,4)
imagesc(CSD_1D);

set(gca, 'yTickLabel', {'DM', 'VL'}, 'yTick', [1 length(diag(squeeze(CSD_est(iS,:,:))))])
set(gca, 'xTickLabel', filt_data.tvec)
colorbar

csd_struct.csd_1d = CSD_1D;
csd_struct.csd_est = CSD_est;
csd_struct.X = X;
csd_struct.Y = Y;
csd_struct.samples = samples;
