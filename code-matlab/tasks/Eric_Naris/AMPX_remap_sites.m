function data_out = AMPX_remap_sites(data, ExpKeys)
%% AMPX_map: rempas the data channels to match the probe layout specified
% in the Expkeys.
%
%  Inputs:
%   - data: [struct] format (AMPX, FT, TSD) will be autodetected
%   - ExpKeys: [struct] contains the mapping
%
%  Outputs:
%   - data_out: [struct] in the same format as the input
%
%
%% determine the data type

if isfield(data, 'trial')
    type = 'FT';
elseif isfield(data, 'channels')
    type = 'AMPX';
elseif isfield(data, 'tvec') && isfield(data, 'data')
    type = 'TSD';
else
    error('Data is not in an acceptable format')
end

Site_map = reshape(ExpKeys.Probe_layout,8,8)

data_out = data;
%% remap the data
switch type
    case 'TSD'
        disp('Data is in TSD format')
        for iChan = length(data.channels):-1:1
            data_out.channel{iChan} = data.channels{ExpKeys.Probe_layout(iChan)};
            data_out.label_remap(iChan) = ExpKeys.Probe_layout(iChan);
            sites_out(iChan) = iChan;
            sites_in(iChan) = ExpKeys.Probe_layout(iChan);
        end
        
        
    case 'FT'
        disp('Data is in FT format')
        for iChan = length(data.label):-1:1
            data_out.trial{1}(iChan,:) = data.trial{1}(ExpKeys.Probe_layout(iChan),:);
            data_out.label_remap{iChan} = ExpKeys.Probe_layout(iChan);
            sites_out(iChan) = iChan;
            sites_in(iChan) = ExpKeys.Probe_layout(iChan);
        end
        
        
    case 'AMPX'
        disp('Data is in AMPX format')
        for iChan = 64:-1:1
            data_out.channels{iChan} = data.channels{ExpKeys.Probe_layout(iChan)};
            data_out.labels_remap(iChan) = ExpKeys.Probe_layout(iChan);
            sites_out(iChan) = iChan;
            sites_in(iChan) = ExpKeys.Probe_layout(iChan);
        end
end

disp('Sites In')
reshape(sites_in, 8,8)
disp('Sites Out')
reshape(sites_out,8,8)