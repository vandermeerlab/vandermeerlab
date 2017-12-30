function S = LoadFrankSpikes(cfg,filename,epoch)
%LOADFRANKSPIKES Loads spike-sorted data provided by Frank and colleages 
% (from CRCNS)
%   S = LoadFrankSpikes(cfg,filename)
%
%   CFG OPTIONS
%     cfg.verbose = 1  (default)
%          1 - Display text in command window
%          0 - Do not display text in command window
%
%   INPUTS
%         cfg - Config struct with fields controlling function behaviour
%    filename - String with extentsion specifying the the spike data you want
%               to load (e.g. 'bonspikes07.mat'); only loads one file at 
%               a time.
%       epoch - Single integer specifying the epoch of interest
%
%   OUTPUTS
%         S - spike-sorted data presented in a TS datatype format
%
% aacarey Dec 2017 initial version
% see also: LoadSpikes

% Declare global variables
global name

cfg_def.verbose = 1;
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);
cfg.epoch = epoch; % keep as a record even though it's not a config option at this time

% Make sure user has specified a file that exists in the directory
if isempty(FindFiles(filename))
    error('File does not exist in directory. Check spelling.');
end

% Get base file name
[~,name,~] = fileparts(filename);

% Don't load the wrong type of data
if ~strcmp(name(4:9),'spikes')
    error('The file is not spike-sorted data')
end



% Talk to me
if cfg.verbose;fprintf('%s: loading %s, epoch %d...\n',mfun,filename,epoch); end

% Load the data
spikes = load(filename,'spikes');

% Get today's spikes
recordingDay = str2double(name(10:11));
recordingDaySpikes = spikes.spikes{1,recordingDay};

% Get all spikes from the epoch
epochSpikes = recordingDaySpikes{1,epoch};

[t,label,TT_num] = getTTSpikes(epochSpikes);

% Make output
S = ts2(t,label);
S.usr.TT_num = TT_num;
S.cfg.SessionID = name;
S = History(S,mfun,cfg);

% Make sure TSD is well formed
if ~CheckTS(S)
    error('Failed to load spike data as a proper TS datatype')
end

if cfg.verbose; fprintf('%s: loaded %d units \n',mfun,length(S.t)); end

%------------HELPERS---------------

    function [t,label,TT_num] = getTTSpikes(epochSpikes)
        nTetrodes = length(epochSpikes);
        nthUnit = 1;
        for iTT = 1:nTetrodes
            for iUnit = 1:length(epochSpikes{1,iTT})
                unitSpikes = epochSpikes{1,iTT}{1,iUnit};
                
                if ~isempty(unitSpikes) && ~isempty(unitSpikes.data)
                    t{1,nthUnit} = unitSpikes.data(:,1);                    
                    label{1,nthUnit} = [name,'-TT',num2str(iTT),'_',num2str(iUnit)];
                    TT_num(nthUnit) = iTT; 
                    nthUnit = nthUnit +1;
                end                
            end
        end
    end
end

