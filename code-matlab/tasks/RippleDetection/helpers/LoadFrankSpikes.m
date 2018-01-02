function S = LoadFrankSpikes(cfg,epoch)
%LOADFRANKSPIKES Loads spike-sorted data provided by Frank and colleages 
% (from CRCNS)
%   S = LoadFrankSpikes(cfg,filename)
%
%   CFG OPTIONS
%     cfg.area = '' (default)
%          'CA1' - Load cells recorded in CA1
%          'CA3' - Load cells recorded in CA3
%
%     cfg.verbose = 1  (default)
%          1 - Display text in command window
%          0 - Do not display text in command window
%
%   INPUTS
%         cfg - Config struct with fields controlling function behaviour
%       epoch - Single integer specifying the epoch of interest
%
%   OUTPUTS
%         S     - spike-sorted data presented in a TS datatype format
%     FIELDS of S:
%          .t   - {1xnCells} cell array containing timestamps of spiking units
%   .usr.TT_num - [1xnCells] double containing the tetrode number the cell
%                 was recorded on
%   .usr.area   - {1xnCells} cell array containing the brain region the
%                 cell was recorded from. (This field exists only if the
%                 animal's cellinfo.mat file is available in the same
%                 directory as the spikes.mat file.)
%
% aacarey Dec 2017 initial version
%   - edit, added usr.area, function automatically looks for file
% see also: LoadSpikes

% Parse config
cfg_def.verbose = 1;
cfg_def.area = '';
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);
cfg.epoch = epoch; % keep as a record even though it's not a config option at this time

% Find spikes file

filename = FindFiles('*spikes*');
if isempty(filename)
    error('spikes file does not exist in directory.');
elseif numel(filename) > 1
    error('More than one spikes file found in directory.')
else
    filename = filename{1}; % pull the thing out of the cell array
end

% Get base file name
[~,name,~] = fileparts(filename);

% Talk to me
if cfg.verbose;fprintf('%s: loading %s, epoch %d...\n',mfun,[name,'.mat'],epoch); end

% Load the data
spikes = load(filename,'spikes');

% Get cell info for usr field 'area'
file = FindFiles(['*',name(1:3),'cellinfo.mat']);
if isempty(file)
    % Two types of warnings in case someone has warning off
    warning('Could not find cellinfo in file directory')
    fprintf('WARNING in %s: Could not find cellinfo.mat. Some functions require usr fields from cellinfo.\n',mfun)
    cellinfo = {};
elseif numel(file) > 1
    error('More than one cellinfo file was found in the current directory')
else
   load(file{1,1}) % var loads as 'cellinfo'
end

% Get today's spikes
recordingDay = str2double(name(10:11));
recordingDaySpikes = spikes.spikes{1,recordingDay};

% Get all spikes from the epoch
epochSpikes = recordingDaySpikes{1,epoch};

nthUnit = 1;
for iTT = 1:length(epochSpikes) % nTetrodes
    for iUnit = 1:length(epochSpikes{1,iTT})
        unitSpikes = epochSpikes{1,iTT}{1,iUnit};
        
        if ~isempty(unitSpikes) && ~isempty(unitSpikes.data)
            t{1,nthUnit} = unitSpikes.data(:,1);
            label{1,nthUnit} = [name,'-TT',num2str(iTT),'_',num2str(iUnit)];
            usr.TT_num(nthUnit) = iTT;
            
            if ~isempty(cellinfo) % put recording area ('CA1' or 'CA3' for example) into usr field
                usr.area{nthUnit} = GetFrankInfo(cellinfo,'area',recordingDay,epoch,iTT,iUnit);
            end
            
            nthUnit = nthUnit +1;
        end
    end
end

% Make output
S = ts2(t,label);
S.usr = usr;

% Select cells recorded from specified brain region
fromArea = '';
if ~isempty(cfg.area)
    fromArea = [' from area ',cfg.area];
    cfg_temp = []; cfg_temp.str = cfg.area; cfg_temp.verbose = 0;
    S = SelectTS(cfg_temp,S,'area');
end

% Keep a record of config history
S.cfg.SessionID = name;
S = History(S,mfun,cfg);

% Make sure TSD is well formed
if ~CheckTS(S)
    error('Failed to load spike data as a proper TS datatype')
end

if cfg.verbose; fprintf('%s: loaded %d units%s \n',mfun,length(S.t),fromArea); end

end

