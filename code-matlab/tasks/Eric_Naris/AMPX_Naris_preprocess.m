function [data, data_ft] = AMPX_Naris_preprocess(cfg_in, file_name, session_type)
%% AMPX_Naris_preprocess: Loads and removees the DC offset, and wide band passes the data for all data.  
%The data has already been preprocessed it will simply load the data.
% Must be in the data folder for this to work. 
%
%   Inputs: 
%        - fname: [string] in the 'R0XX-yyyy-mm-dd' format
%        - session_type: [string] 'pre', 'task', or 'post'. used to
%        determine which data files to use
%
%   Outputs:
%        - data: [struct] in AMPX format contains
%                      - hdr
%                      - channels
%                      - labels
%                      - tvec
%                      - dc_remove: 
%
%
%% configs
cfg_def.wideband = [1 500];
cfg_def.buttord = 10;
cfg_def.dec_fac = 10;
cfg = ProcessConfig(cfg_def, cfg_in);

%%
if exist([session_type '_data_preprocess_filt.mat'],'file')==0 % check to see if the data has already been preprocessed. If so, load it. 
    if exist([ session_type '_data_preprocess.mat'],'file')~=0 % if the data has been loaded but not wideband filtered
        disp('Data file found')
        load([session_type '_data_preprocess.mat'])
        disp('Data file loaded')
    else
        disp('No data file found..')
        disp('Loading Data');
        if strcmp(session_type, 'task') == 1;
            channels_to_load = 1:74;
        else
            channels_to_load = 1:64;
        end
        data = AMPX_loadData([file_name '.dat'],channels_to_load,cfg.dec_fac); % note decimation factor of 20
        disp('Data loading complete');
        save([session_type '_data_preprocess.mat'], 'data', '-v7.3')
    end
    session_name = [data.hdr.Filename(5:8) '_' data.hdr.Filename(10:13) '_' data.hdr.Filename(15:16) '_' data.hdr.Filename(18:19)];
    
    %% DC remove (only if this has not previously been done.
    if isfield(data,'data_ft') ==0
        for iChan = size(data.channels,2):-1:1
            data_temp{1,iChan} = data.channels{1,iChan} - mean(data.channels{1,iChan});
        end
        data.channels = data_temp;
        data.dc_remove = 'yes';
        clear data_temp;

        %% Bandpass the signal
        [data.data_ft] = AMPX_filter(data, cfg.wideband(1), cfg.wideband(2), cfg.buttord);
        
        save([session_type '_data_preprocess_filt.mat'], 'data', '-v7.3')
    end
else
    disp(['Data has previously been preprocessed, loading data'])
    load([session_type '_data_preprocess_filt.mat'])
    session_name = [data.hdr.Filename(5:8) '_' data.hdr.Filename(10:13) '_' data.hdr.Filename(15:16) '_' data.hdr.Filename(18:19)];
end

disp(['Data loaded: ' data.hdr.Filename])
%% pull the data type out of the data struct
data_ft = data.data_ft;
data = rmfield(data, 'data_ft');