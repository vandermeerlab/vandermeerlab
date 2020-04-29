%%
restoredefaultpath;
addpath('D:\My_Documents\GitHub\fieldtrip'); ft_defaults
addpath('D:\My_Documents\GitHub\striatal-spike-rhythms\shared\io\ft');
addpath('D:\My_Documents\GitHub\striatal-spike-rhythms\shared\io\neuralynx');

%%
cd('D:\data\adrlab\R117\R117-2007-06-10');

run(ls('*keys.m'));

fc = ExpKeys.goodGamma_vStr;
data_all = ft_read_neuralynx_interp(fc);

%% restrict to on-task data only

t0 = double((data.hdr.FirstTimeStamp ./ data.hdr.TimeStampPerSample)) / data.hdr.Fs; % time of first sample in ft data
t1 = ExpKeys.TimeOnTrack - t0; t2 = ExpKeys.TimeOffTrack - t0;
cfg = []; cfg.latency = [t1 t2];
data_all = ft_selectdata(cfg, data_all);

%% STS (basic)
% cfg              = [];
% cfg.method       = 'mtmfft';
% cfg.foilim       = [1 100]; % cfg.timwin determines spacing
% cfg.timwin       = [-0.5 0.5];
% cfg.taper        = 'hanning';
% cfg.spikechannel = data_all.label(2);
% cfg.channel      = data_all.label(1);
% stsFFT           = ft_spiketriggeredspectrum(cfg, data_all);

%% STS with spike times (not binarized)
fc = dir('*.t'); fc = {fc.name};
for iF = 5%1:length(fc)
    spike = ft_read_spike(fc{iF}); % needs fixed read_mclust_t.m (edited D:\My_Documents\GitHub\fieldtrip\fileio\private\read_mclust_t.m to uint32)
    
    cfg           = [];
    cfg.hdr       = data_all.hdr; % contains information for conversion of samples to timestamps
    cfg.trlunit   = 'samples';
    cfg.trl       = [data_all.sampleinfo(1) data_all.sampleinfo(2) 0]; % now in samples
    spike_trl     = ft_spike_maketrials(cfg,spike);
    
    % only proceed if more than some number of spikes, spike not the same as LFP channel etc..
    
    
    % spike triggered spectrum (nice convolution method)
    cfg              = [];
    cfg.method       = 'mtmconvol';
    cfg.foi          = 1:100;
    cfg.t_ftimwin    = 5./cfg.foi; % 5 cycles per frequency
    cfg.taper        = 'hanning';
    cfg.channel      = data_all.label(1);
    stsConvol        = ft_spiketriggeredspectrum(cfg, data_all, spike_trl);
    
    % compute the statistics on the phases
    cfg               = [];
    cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency
    cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
    cfg.timwin        = 'all'; % compute over all available spikes in the window
    %statSts           = ft_spiketriggeredspectrum_stat(cfg, stsFFT);
    statSts           = ft_spiketriggeredspectrum_stat(cfg, stsConvol);
    
    % plot the results
    figure
    plot(statSts.freq, statSts.ppc0')
    xlabel('frequency')
    ylabel('PPC')
    
end
