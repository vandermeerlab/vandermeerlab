figure
hold on
itrial = 1;
boost = 200;
loop = 1;
iChan = diag(reshape(ExpKeys.Probe_layout,8,8));
c_ord = linspecer(length(iChan));
for ii = 1:length(iChan)
    disp(iChan(ii))
    plot(data.data_ft.time{1}(4000:6000), data.data_ft.trial{1}(iChan(ii),4000:6000) + (boost *ii), 'color', c_ord(ii,:))
    text(2.1, (boost *ii), num2str(iChan(ii)))
end


%%
figure
hold on
itrial = 1;
boost = 200;
loop = 1;
iChan = diag(reshape(ExpKeys.Probe_layout,8,8));
c_ord = linspecer(length(iChan));
for ii = 1:length(iChan)
    disp(iChan(ii))
    plot(fake_AMPX.tvec, fake_AMPX.channels{iChan(ii)} + (boost *ii), 'color', c_ord(ii,:))
    text(2.1, (boost *ii), num2str(iChan(ii)))
end


%%
data_tsd = AMPX_to_tsd(data);
csc = data_tsd;
csc.data = csc.data(ExpKeys.DetectChan,:);
csc.detect_chan = ExpKeys.DetectChan;
    cfg_plot.display = 'tsd'; % tsd, iv
    PlotTSDfromIV(cfg_plot,rand_iv.low,csc);