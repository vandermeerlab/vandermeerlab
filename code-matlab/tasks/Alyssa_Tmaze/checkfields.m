function [tf,missingFields] = checkfields(cfg_in,checkStruct,reqFields)
%CHECKFIELDS Verify the existence of specified fields in a struct.
%   [TF,missingFields] = checkfields(cfg_in,checkStruct,reqFields)
%
%   
%   INPUTS
%             cfg - config with optional fields (see CONFIGS for more
%             details).
%
%     checkStruct - struct whose fields you want to check, ex: ExpKeys
%                 - string containing wildcard input (globfn of FindFiles); 
%                   it's important to have the correct extension 
%                            ex: '*keys.m' or '*metadata.mat'
%                   NOTE: if you are running a script (ext is .m) this 
%                   function assumes it's ExpKeys!
%                   NOTE: a .mat file must contain only one struct for this
%                   function to work.
%     reqFields   - cell array of strings containing field names that are 
%                   required in checkStruct, 
%                      ex: {'TimeOnTrack','TimeOffTrack','nTrials'}
%                      MetaFields = {'coord','taskvars'};
%
%   CONFIGS
%
%     cfg.fd      - string specifying the directory to search in; default
%                   is the current directory (pwd)
%     cfg.verbose - 0 default (do not display missing fields); 1 (display
%                   missing fields). If cfg.verbose = 1 and no fields are
%                   missing, nothing is displayed. 
%
%   OUTPUTS
%
%     tf          - returns 0 if no fields are missing, and 1 if at least
%                   one field is missing. Can be used as a flag. 
%   missingFields - 1 x nMissing cell array containing the names of missing
%                   fields.
%
%
% A.Carey Feb 2015.

%% Parse cfg parameters

cfg_def.verbose = 0; % display the missing fields
cfg = ProcessConfig2(cfg_def,cfg_in);

%%

if isfield(cfg,'fd')
    originalFolder = pwd;
    cd(cfg.fd)
end

if isa(checkStruct,'struct')
        missingFields = reqFields(arrayfun(@(x) ~isfield(checkStruct,x),reqFields));

        if ~isempty(missingFields)
            tf = 1;
        else
            tf = 0;
        end
        
        if cfg.verbose && tf
            [~,sessionName,~] = fileparts(pwd);
            structName = inputname(2);
            disp([sessionName,' is missing ',structName,' field(s):'])
            disp(missingFields)
        end
                      
elseif isa(checkStruct,'char')
    fn = FindFiles(checkStruct);

    if isempty(fn)
        error('No files matching the wildcard input were found.')
    else [~,fileName,ext] = fileparts(fn{1});
    end

    % bring struct into caller workspace

    if strcmp(ext,'.m')
        % it's a script, most likely ExpKeys
        run(fileName);
        checkStruct = ExpKeys;
        structName = 'ExpKeys';
    else
        % it's a mat file
        wubbaffet = load(fileName);
        wubbawubba = fieldnames(wubbaffet);
        checkStruct = getfield(wubbaffet,wubbawubba{1});
        % getting seriously tired of naming variables 
        structName = wubbawubba{1};
    end

    missingFields = reqFields(arrayfun(@(x) ~isfield(checkStruct,x),reqFields));

    if ~isempty(missingFields)
        tf = 1;
    else
        tf = 0;
    end
    
    if cfg.verbose && tf 
        [~,sessionName,~] = fileparts(pwd);
        disp([sessionName,' is missing ',structName,' field(s):'])
        disp(missingFields)
    end

else
    cd(originalFolder)
    error('checkStruct must be a struct or a string that identifies a file containing a struct.')
end

if isfield(cfg,'fd')
    cd(originalFolder)
end

end
