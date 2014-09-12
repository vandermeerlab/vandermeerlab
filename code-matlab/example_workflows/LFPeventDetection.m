%% some hippocampus data
cd('C:\data\R042-2013-08-18_recording');

%% load a CSC (CA1 layer I think)
cfg = []; cfg.fc = {'CSC17.Ncs'};

lfp = LoadCSC(cfg);

%% filter in low gamma band
cfg = [];
cfg.f = [45 55];
cfg.display_filter = 0;

lg = FilterLFP(cfg,lfp);

%% obtain power and z-score it
lgp = LFPpower([],lg);
lgp_z = zscore_tsd(lgp);

%% detect events
cfg = [];
cfg.method = 'raw';
cfg.threshold = 3;
cfg.dcn =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.minlen = 0.05; % minimum interval length

lgp_evt = TSDtoIV(cfg,lgp_z);

%% to each event, add a field with the max z-scored power (for later selection)
cfg = [];
cfg.method = 'max'; % 'min', 'mean'
cfg.label = 'maxlgp'; % what to call this in iv, i.e. usr.label

lgp_evt = AddTSDtoIV(cfg,lgp_evt,lgp_z);

%% select only those events of >5 z-scored power
cfg = [];
cfg.dcn = '>';
cfg.threshold = 5;

lgp_evt = SelectIV(cfg,lgp_evt);

%% plot events in highlighted on top of full lfp
PlotTSDfromIV([],lgp_evt,lfp);

%% ..or the events alone (fixed 200ms window centered at event time)
close all;

cfg = [];
cfg.display = 'iv';
cfg.mode = 'center';
cfg.fgcol = 'k';

PlotTSDfromIV(cfg,lgp_evt,lfp);
%% ..hold on (highlight edges of event on top of previous plot)
cfg = [];
cfg.display = 'iv';
cfg.fgcol = 'r';

PlotTSDfromIV(cfg,lgp_evt,lfp);