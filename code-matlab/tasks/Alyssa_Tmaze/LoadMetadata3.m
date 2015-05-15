function LoadMetadata3(cfg_in)
%LOADMETADATA3 Load metadata struct into caller workspace
%   Use this if you want to load metadata with a different suffix
%        ex: RXXX-20XX-XX-XX-metadata_test.mat, where '_test' is the suffix
%
% The metadata file must contain only one variable called 'metadata' for this
% to work as expected. The variable metadata may contain any number of
% additional fields. The variable can be loaded with a name corresponding
% to its filename if you specify cfg.rename = 1
%
%   LoadMetadata3(cfg_in)
%         cfg_def.suffix = ''; % the suffix for the metadata file you want 
%                 to load (this does not include the extension)
%                 ex: '_test' or '-14-May-2015' thus it will search for
%                 '*metadata_test.mat' or '*metadata-14-May-2015.mat'
%         cfg_def.rename = 0; % don't rename the metadata variable to match 
%                 the filename; 1 rename the variable
%                 ex: loads variable with the name metadata_test or 
%                     metadata_15_May_2015
%
% ACarey May 2015 

%%
cfg_def.suffix = '';
cfg_def.ext = '';
cfg_def.rename = 0; % don't rename the metadata variable to match the filename
cfg = ProcessConfig2(cfg_def,cfg_in);

%[~,sessionID,~] = fileparts(pwd); 
fileID = ['metadata',cfg.suffix];
%fn = FindFiles([sessionID,'-',fileID,'.mat']);
fn = FindFiles(['*',fileID,'.mat']);
if isempty(fn)
    disp(['LoadMetadata3: No files matching',' ''*',fileID,'.mat'' ','were found in', [' ',pwd]])
elseif length(fn) > 1
    disp(['LoadMetadata3: More than one file matching',' ''*metadata.mat'' ','was found in', [' ',pwd]])
    disp('metadata not loaded')
else
    load(fn{1})
    
    if cfg.rename
        % a variable name cannot contain the character '-', so replace with '_'
        for iChar = 1:length(fileID)
            if strcmp(fileID(iChar),'-')
                fileID(iChar) = '_';
            end
        end
        assignin('caller',fileID,metadata)
    else
        assignin('caller','metadata',metadata)
    end
end

end

