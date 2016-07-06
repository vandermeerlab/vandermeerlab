function ExpKeys_remap = AMPX_remap_ExpKeys(ExpKeys_in, data_remap_AMPX)
%% ExpKeys_remap: remaps the channel labels for the "probe_layout" field      

%          Inputs: 
%           - ExpKeys [struct]
%          Outputs: 
%           - ExpKeys_remap [struct]
% 
% EC - 2016-05-18

%% 
ExpKeys_remap = ExpKeys_in;

% convert the probe layout
ExpKeys_remap.Probe_layout = reshape(data_remap_AMPX.labels(1:64),8,8);

ExpKeys_remap.DetectChan = find(ExpKeys_in.Probe_layout == Naris_BestChan_remap(ExpKeys_in, 'location','vl'));
% convert the bad channels as well. 
for ii = 1:length(ExpKeys_in.BadChannels)
    ExpKeys_remap.BadChannels(ii) = find(reshape(ExpKeys_in.Probe_layout,8,8) == ExpKeys_in.BadChannels(ii));
end

disp('ExpKeys Channels remapped')
ExpKeys_remap.remapped = 'yes';
end