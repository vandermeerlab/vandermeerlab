function tsd_out = AMPX2tsd(data)


for ichan = length(data.channels):-1:1
    data_temp(ichan,:) = data.channels{ichan};
end

tsd_out = tsd(data.tvec, data_temp, data.labels);
for ichan = length(data.channels):-1:1
tsd_out.cfg.hdr{ichan}.SamplingFrequency = data.hdr.Fs;
end