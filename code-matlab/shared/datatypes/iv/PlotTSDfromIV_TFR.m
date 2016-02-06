function PlotTSDfromIV_TFR(cfg_in,iv_in,ft_in)
% function PlotTSDfromIV_TFR(cfg,iv_in,ft_in)
%
% plot ft_in data spectrograms and raw data as defined by iv_in
%
% INPUTS:
%
% ft_in: input fieldtrip data to plot from
% iv: intervals to plot data for
%
% CFG OPTIONS:
%
% cfg.target = []; % if ft_in has multiple data dimensions (labels)
% cfg.twin = [-0.5 0.5];
% cfg.clim = [lower_lim upper_lim] used to make the caxis in the imagesc plots even. default is [0 500]
% ...and a bunch more
%
% MvdM 2014-11-12

cfg_def              = [];
cfg_def.twin         = [-0.5 0.5];
cfg_def.dt           = 0.01;
cfg_def.method       = 'mtmconvol';
cfg_def.taper        = 'hanning';
cfg_def.foi          = 1:1:100; % frequencies of interest
cfg_def.subplotdim   = [4 5];
cfg_def.clim         = [0 500]; % sets the caxis for the imagesc plots. default [0 500]

cfg = ProcessConfig2(cfg_def,cfg_in); % should take whatever is in cfg_in and put it into cfg!

cfg.t_ftimwin    = 15*(1./cfg.foi);
cfg.toi          = cfg.twin(1):cfg.dt:cfg.twin(2);

% check if conditions are in place
nTrials = length(ft_in.trial);
if nTrials > 1
   error('Multiple trials not yet implemented.');
   % MvdM: ..and may not even make sense
end

nData = size(ft_in.trial{1},1);
if nData > 1
    if ~isempty(cfg.target)
        temp_cfg = [];
        temp_cfg.channel = cfg.target;
        ft_in = ft_selectdata(temp_cfg,ft_in);
    else
       error('Multiple data dimensions exist but no label is specified.'); 
    end
end

% create trials from iv
trl_cfg = [];
trl_cfg.t = IVcenters(iv_in);
trl_cfg.mode = 'ft';
trl_cfg.hdr = ft_in.hdr;

trl_cfg.twin = cfg.twin + [-2 2]; % padding for spectrogram
trl_cfg.hdr.FirstTimeStamp = ft_in.time{1}(1); 
trl_cfg.tvec = ft_in.time{1}; % For neuralynx data conversion

trl = ft_maketrl(trl_cfg);

temp_cfg = []; temp_cfg.trl = trl;
ft_in = ft_redefinetrial(temp_cfg,ft_in);

% specgram
cfg.output       = 'pow';
cfg.keeptrials   = 'yes'; % need this for stats later

TFR = ft_freqanalysis(cfg, ft_in);

% plot
temp_cfg = [];
temp_cfg.latency = cfg.twin;
TFR = ft_selectdata(temp_cfg,TFR);
ft_in = ft_selectdata(temp_cfg,ft_in);

nppf = prod(cfg.subplotdim);
for iI = 1:length(ft_in.trial)
    
    figno = ceil(iI./nppf);
    plotno = mod(iI-1,nppf) + 1;
    maximize;
    
    figure(figno);
    subtightplot(cfg.subplotdim(1),cfg.subplotdim(2),plotno);
    
    imagesc(TFR.time,TFR.freq,sq(TFR.powspctrm(iI,1,:,:))); axis xy
    caxis(cfg.clim)
    hold on;
    
    temp_data = rescale(ft_in.trial{iI},cfg.foi(1),cfg.foi(end));
    % note, may need alternative option that doesn't rescale
    plot(ft_in.time{iI},temp_data,'LineWidth',1,'Color',[1 1 1]);

    axis off; axis tight; 
end