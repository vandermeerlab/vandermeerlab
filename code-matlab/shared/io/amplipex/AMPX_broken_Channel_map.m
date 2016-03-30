function AMPX_broken_Channel_map(ExpKeys)
for iChan = 1:64
    if ismember(ExpKeys.Probe_layout(iChan), ExpKeys.BadChannels)
        Map(iChan) = NaN;
    else
        Map(iChan) = ExpKeys.Probe_layout(iChan);
    end
end

reshape(Map,8,8)