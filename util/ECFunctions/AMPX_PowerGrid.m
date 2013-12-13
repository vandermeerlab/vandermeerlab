
function [data] = AMPX_PowerGrid(data,varargin)
%AMPX_PowerGrid will generate an 8x8 matrix of the power for each frequency range of interest as well as heat maps for each frequency
%
%
%INPUTS:
%
% data: [1x1 struct] - this is the output of AMPX_loadData containing the
% hdr, channels, tvec, and labels.  This must have been sorted based on the
% physical location of each recording site.  This can be done using the
% AMPX_Sort_Channels function and replacing the data.channels with the
% sorted_data cells given by AMPX_sort_Channels.
%
%
% OUTPUTS:
%
% data.power: [1 x 1 struct] Adds an additional field to the data structure
% that contains the average power across the delta,theta,beta,gamma_low and
% gamma_high bands.  It will also plot the heat maps for each frequency
% band across all the recording sites.

%% Initialize variables
Window = 2; %in seconds should be able to adapt to the decimation
windowSize = 4096;
% Overlap
% NFFt must be power of two
extract_varargin;

%% Prepare the smapling window values.
nFFT_vec = 1:32;
nFFT_vec = 2.^nFFT_vec;
nFFT_ind = nFFT_vec > windowSize;
%%  Extract the PSD and then create a master plot for all the channels
PSD_all = figure(1);
data.PSD = cell(1,64); data.Freq = cell(1,64);
for iChan = 1:size(data.channels,2)
    subplot(8,8,iChan)
    [Pxx,F] = pwelch(data.channels{iChan},hanning(4096),3072,16384,(data.hdr.Fs));
    data.PSD{1,iChan} = 10*log10(abs(Pxx)); data.Freq{1,iChan} = F;
    plot(F,10*log10(abs(Pxx)),'b'); xlim([0 200]);
    xlabel('Frequency (Hz)'); ylabel('power'); title('PSD');
    title(data.labels(iChan));
end

% add other window sizes like 4 s, 8s

%% Extract the mean power for each band of interest and place them in a grid for each site.

for iChan = 1:size(data.channels,2)
    data.power.gamma_low(1,iChan) = nanmean(data.PSD{1,iChan}(data.Freq{1,iChan}>=45 & data.Freq{1,iChan}<=55));  % finds the low Gamma power
    data.power.gamma_high(1,iChan) = nanmean(data.PSD{1,iChan}(data.Freq{1,iChan}>=70 & data.Freq{1,iChan}<=100));  % finds the high Gamma power
    data.power.theta(1,iChan) = nanmean(data.PSD{1,iChan}(data.Freq{1,iChan}>=7 & data.Freq{1,iChan}<=10));  % finds the Theta power
    data.power.beta(1,iChan) = nanmean(data.PSD{1,iChan}(data.Freq{1,iChan}>=15 & data.Freq{1,iChan}<=25));  % finds the Beta power
    data.power.delta(1,iChan) = nanmean(data.PSD{1,iChan}(data.Freq{1,iChan}>=1.5 & data.Freq{1,iChan}<=3));  % finds the delta power
end

% also grab a psd mean of the entire PSd and the 1-200 HZ band.  This is to
% ensure that this is not some effect of the entire spectra.

% Set the values of known broken channels to NaN
data.power.gamma_low(data.bad_channels)=NaN;
data.power.gamma_high(data.bad_channels)=NaN;
data.power.beta(data.bad_channels)=NaN;
data.power.theta(data.bad_channels)=NaN;
data.power.delta(data.bad_channels)=NaN;

% form a matrix with each shank as a column.
data.power.gamma_low = reshape(data.power.gamma_low,8,8);
data.power.gamma_high = reshape(data.power.gamma_high,8,8);
data.power.beta = reshape(data.power.beta,8,8);
data.power.theta = reshape(data.power.theta,8,8);
data.power.delta = reshape(data.power.delta,8,8);


%% make the plots.
% Delta
Power_delta = figure(2);
imagesc(data.power.delta)
grid on
axh = gca;
set(axh,'XTick',[0.5:1:8.5])
set(axh,'yTick',[0.5:1:8.5])
xlabel('Ventral'); ylabel('Medial');
title(['Delta: R036-' data.hdr.date '-' data.hdr.Filename(end-6:end-4)])
y = 0;  z = 0;
for ii = 1:8
    y = y+1;
    x = 0;
    for jj = 1:8
        x = x+1; z = z+1;
        th = text(ii-.2,jj,num2str(data.labels(z))); set(th,'Color',[1 1 1],'FontSize',24)
        % disp(z)
    end
end
saveas(Power_delta,['D:\DATA\R036\PSD_plots\' data.hdr.Filename(5:end-4) 'Delta'],'png')


%% Theta
Power_theta = figure(3);
imagesc(data.power.theta)
grid on
axh = gca;
set(axh,'XTick',[0.5:1:8.5])
set(axh,'yTick',[0.5:1:8.5])
xlabel('Ventral'); ylabel('Medial');
title(['Theta: R036-' data.hdr.date '-' data.hdr.Filename(end-6:end-4)])
y = 0;  z = 0;
for ii = 1:8
    y = y+1;
    x = 0;
    for jj = 1:8
        x = x+1; z = z+1;
        th = text(ii-.2,jj,num2str(data.labels(z))); set(th,'Color',[1 1 1],'FontSize',24)
        % disp(z)
    end
end
saveas(Power_theta,['D:\DATA\R036\PSD_plots\' data.hdr.Filename(5:end-4) 'Theta'],'png')

%% Beta
Power_beta = figure(4);
imagesc(data.power.beta)
grid on
axh = gca;
set(axh,'XTick',[0.5:1:8.5])
set(axh,'yTick',[0.5:1:8.5])
xlabel('Ventral'); ylabel('Medial');
title(['Beta: R036-' data.hdr.date '-' data.hdr.Filename(end-6:end-4)])
y = 0;  z = 0;
for ii = 1:8
    y = y+1;
    x = 0;
    for jj = 1:8
        x = x+1; z = z+1;
        th = text(ii-.2,jj,num2str(data.labels(z))); set(th,'Color',[1 1 1],'FontSize',24)
        % disp(z)
    end
end
saveas(Power_beta,['D:\DATA\R036\PSD_plots\' data.hdr.Filename(5:end-4) 'Beta'],'png')

%% Low Gamma
Power_gamma_low = figure(5);
imagesc(data.power.gamma_low)
grid on
axh = gca;
set(axh,'XTick',[0.5:1:8.5])
set(axh,'yTick',[0.5:1:8.5])
xlabel('Ventral'); ylabel('Medial');
title(['Low Gamma: R036-' data.hdr.date '-' data.hdr.Filename(end-6:end-4)])
y = 0;  z = 0;
for ii = 1:8
    y = y+1;
    x = 0;
    for jj = 1:8
        x = x+1; z = z+1;
        th = text(ii-.2,jj,num2str(data.labels(z))); set(th,'Color',[1 1 1],'FontSize',24)
        % disp(z)
    end
end
saveas(Power_gamma_low,['D:\DATA\R036\PSD_plots\' data.hdr.Filename(5:end-4) 'Low_Gamma'],'png')

%% High Gamma
Power_gamma_high = figure(6);
imagesc(data.power.gamma_high)
grid on
axh = gca;
set(axh,'XTick',[0.5:1:8.5])
set(axh,'yTick',[0.5:1:8.5])
xlabel('Ventral'); ylabel('Medial');
title(['High Gamma: R036-' data.hdr.date '-' data.hdr.Filename(end-6:end-4)])
y = 0;  z = 0;
for ii = 1:8
    y = y+1;
    x = 0;
    for jj = 1:8
        x = x+1; z = z+1;
        th= text(ii-.2,jj,num2str(data.labels(z))); set(th,'Color',[1 1 1],'FontSize',24)
        % disp(z)
    end
end
saveas(Power_gamma_high,['D:\DATA\R036\PSD_plots\' data.hdr.Filename(5:end-4) 'High_Gamma'],'png')

end