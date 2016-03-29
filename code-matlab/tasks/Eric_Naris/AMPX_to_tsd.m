function data_tsd = AMPX_to_tsd(data)
%% AMPX_to_tsd: converts the AMPX data format to the tsd format
%
%     Inputs:
%       - data [struct]: AMPX format
%
%     Output:
%       - data_tsd [struct]: tsd format 
%
%% Convert from AMPX to TSD
    for iChan = length(data.channels):-1:1
        data.tsd_channels(iChan,:) = data.channels{iChan};
    end       
    data_tsd = tsd(data.tvec',data.tsd_channels);
    data_tsd.cfg.hdr{1}.SamplingFrequency = data.hdr.Fs;
    data_tsd.label = data.labels;
    if isfield(data, 'labels_remap')
            data_tsd.label_remap = data.labels_remap;
    end