function iv = AddTSDtoIV(cfg_in,iv,tsd_in)
% function iv_out = AddTSDtoIV(cfg_in,iv_in,tsd_in)
%
% add usr field to iv based on tsd
%
% INPUTS:
%
% iv1: interval data to be selected from based on..
% iv2: interval data to be used to select from iv1
%
% CFG OPTIONS:
% cfg.method = 'max'; 'min', 'mean'
% cfg.target = []; % which label in tsd to use
% cfg.label = []; % what to call this in iv, i.e. usr.label
%
% OUTPUTS:
%
% iv_out: output interval data
%
% MvdM 2014-06-24

cfg_def.method = 'max';

cfg = ProcessConfig(cfg_def,cfg_in); % should take whatever is in cfg_in and put it into cfg!
mfun = mfilename;

% check for well formed input
nData = size(tsd_in.data,1);
if nData == 1
    data_in = tsd_in.data;
else
    if isempty(cfg.target)
        error('Multiple data dimensions present, must specify cfg.target.');
    else
        data_in = getd(tsd_in,cfg.target);
    end
end

% find indices for iv
tstart_idx = nearest_idx3(iv.tstart,tsd_in.tvec);
tend_idx = nearest_idx3(iv.tend,tsd_in.tvec);

% collect data
data_temp = nan(size(tstart_idx));

switch cfg.method
    case 'max'
        for iI = 1:length(data_temp)
            data_temp(iI) = max(data_in(tstart_idx(iI):tend_idx(iI)));
        end
    case 'mean'
        for iI = 1:length(data_temp)
            data_temp(iI) = mean(data_in(tstart_idx(iI):tend_idx(iI)));
        end
end

iv.usr.(cfg.label) = data_temp(:); % ensure column shape

% housekeeping
iv.cfg.history.mfun = cat(1,iv.cfg.history.mfun,mfun);
iv.cfg.history.cfg = cat(1,iv.cfg.history.cfg,{cfg});