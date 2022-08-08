function AMPX_imagesc_out(data, nplots)
%% shows a series of imagesc's for the input data :
% 
%
% example: AMPX_imagesc(data_out.random_low_gamma.power.power_distrib, 16)
% will give an 4x4 grid of the imagesc
%% 
figure
hold on
for ii = 1:nextpow2(nplots)^2
    subplot(nextpow2(nplots), nextpow2(nplots), ii)
    nan_imagesc_ec(data{ii});
%     colorbar('AxisLocation','in')
    xlabel([num2str(ii) '    Range: ' num2str(round(min(min(data{ii})))) '-' num2str(round(max(max(data{ii}))))])
end
% tightfig