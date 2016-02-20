function pos_tsd = LoadPos(cfg_in)
% pos_tsd = LoadPos(cfg)
%
% Wrapper function for loading Neuralynx .nvt file
%
% input cfg fields:
%
% cfg.fn: filename to load
%   if no filename is specified, loads *.nvt file in current dir
% cfg.tsflag = 'sec';
% cfg.removeZeros = 1;
% cfg.convFact = [xConvFact yConvFact]; % output from PosCon() function,
%                converts position data units from pixels to centimeters
% cfg.verbose = 1; 1 display command window text, 0 don't
%
% output:
%
% pos_tsd: position data tsd struct
%
% MvdM 2014-06-17
% youkitan 2014-11-05
% A.Carey 2015-02-11 (added convFact)

cfg_def.fn = {};
cfg_def.removeZeros = 1;
cfg_def.tsflag = 'sec';
cfg_def.convFact = [];
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % note: this converts cfg fields into workspace variables!

if isempty(cfg.fn)
    
    cfg.fn = FindFile('*.nvt');
    
end

% load data
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT(cfg.fn, [1 1 1 1 1 1], 1, 1, [] );

switch cfg.tsflag
    case 'sec'
        Timestamps = Timestamps*10^-6;
    case 'timestamp'
        Timestamps = Timestamps;
end

keep_idx = find(X ~= 0 & Y ~= 0);
if cfg.removeZeros

    X = X(keep_idx);
    Y = Y(keep_idx);
    Timestamps = Timestamps(keep_idx);
    
end
if cfg.verbose; fprintf('LoadPos.m: %.2f%% of samples tracked.\n',(length(keep_idx)./length(X)).*100); end
    
    
% fill tsd
pos_tsd = tsd;

pos_tsd.tvec = Timestamps;

pos_tsd.data(1,:) = X;
pos_tsd.label{1} = 'x';

pos_tsd.data(2,:) = Y;
pos_tsd.label{2} = 'y';

% Convert data from pixels to cm
if ~isempty(cfg.convFact) 
    pos_tsd.data(1,:) = pos_tsd.data(1,:)./cfg.convFact(1); 
    pos_tsd.data(2,:) = pos_tsd.data(2,:)./cfg.convFact(2); 
end

% check if ExpKeys available
keys_f = FindFiles('*keys.m');
if ~isempty(keys_f)
    run(keys_f{1});
    pos_tsd.cfg.ExpKeys = ExpKeys;
end

% add sessionID
[~,pos_tsd.cfg.SessionID,~] = fileparts(pwd);

% housekeeping
pos_tsd.cfg.history.mfun = cat(1,pos_tsd.cfg.history.mfun,mfun);
pos_tsd.cfg.history.cfg = cat(1,pos_tsd.cfg.history.cfg,{cfg});