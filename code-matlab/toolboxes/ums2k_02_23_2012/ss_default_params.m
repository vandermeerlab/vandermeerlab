
function spikes = ss_default_params(Fs, varargin )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_default_params - initializes a spikes object with the default
%                     parameters
%               
% Usage:
%       spikes = ss_default_params(Fs, varargin )
%
% Description:  
%   This function creates an empty spikes object containing a params struc.
% There are 2 major categories.  All fields in SPIKES.PARAMS are related
% to the algorithmic process of spike detection and sorting.  In 
% SPIKES.PARAMS.DISPLAY there are parameters that largely only affect the
% way in which data is plotted.
%
%  The only required input is Fs, the data sampling rate in Hertz.  Any
%  other parameters can be changed by putting in a parameter/value pair 
%  such as.  spikes = ss_default_params(36000, 'thresh', 4.0);
%

    %%
    % ALGORITHMIC PARAMETERS
    %
  
    % spike detection parameters
    spikes.params.Fs = Fs;                 % Hz, sampling rate of spike data
    spikes.params.detect_method = 'auto';  % 'auto' = threshold calculated from background noise, 'manual' = user defined threshold
    spikes.params.thresh = 3.9;            % for 'auto', set number of standard deviations above background noise
 %  spikes.params.thresh =    [-5.3897 -5.0999 -5.1816 -5.5274]; % for 'manual', set specific thresholds for each channel
    spikes.params.window_size  = 1.5;       % ms, width of a spike
    spikes.params.shadow = 0.75;            % ms, enforced dead region after each spike
    spikes.params.cross_time = 0.6;         % ms, alignment point for peak of waveform
 
    % sorting parameters
    spikes.params.refractory_period = 2.5;  % ms, refractory period (for calculation refractory period violations)
    spikes.params.max_jitter = 0.6;         % ms, width of window used to detect peak after threshold crossing
    spikes.params.agg_cutoff = .05;         %  higher = less aggregation, lower = more aggregation
    spikes.params.kmeans_clustersize = 500; %  target size for miniclusters

    %%
    % DISPLAY PARAMETERS
    %
          
    % plot_waveforms
    spikes.params.display.default_waveformmode = 3;  % 1 = show all waveforms, 2 = show 95% waveform bands, 3 = show 2d histogram
    spikes.params.display.time_scalebar = 1;         % ms, length of scalebar in milliseconds on waveform plots  
    spikes.params.display.cmap          = hot(64);      % default color map
    
    % plot_features
    spikes.params.display.xchoice = 'PC';         % what feature to show by default (x-axis), see "plot_features.m"
    spikes.params.display.xparam = 1;             % feature parameter (x-axis)
    spikes.params.display.ychoice = 'PC';         % what feature to show by default (y-axis)
    spikes.params.display.yparam = 2;             % feature parameter (y-axis)
    spikes.params.display.show_outliers = 1;      % 0 = hide outliers, 1 = show outliers

    % spike-train correlations parameters
    spikes.params.display.show_isi = 0;                 % 0 = show autocorrelation by default, 1 = show ISI histogram
    spikes.params.display.max_autocorr_to_display = .1; % s, maximum time lag in autocorrelation/cross-correlation plots
    spikes.params.display.max_isi_to_display = .025;    % s, maximum time lag to show in ISI histograms,
    spikes.params.display.correlations_bin_size = 2;    % ms, bin width used for spike time correlations
    spikes.params.display.isi_bin_size = .5;            % ms, bin width used for ISI histograms
    spikes.params.display.default_xcorr_mode = 1;       % 0 = show autocorrelation of merged spike trains by default, 1 = show x-correlation 
    spikes.params.display.trial_spacing = .5;           % s, time padding between trials
    
    % plot_stability parameters
    spikes.params.display.stability_bin_size = 10;      % s, bin used for estimating firing rate
    spikes.params.display.max_scatter = 1000;           % maximum number of waveforms to show in amplitude scatter plot

    % outlier tool
    spikes.params.display.default_outlier_method = 1; % 1 = estimate covariance from cluster, 2 = estimate covariance from background noise        
    
    % cluster labeling in splitmerge tool, the first option is the default
    spikes.params.display.label_categories = {'in process', 'good unit', 'multi-unit', 'garbage', 'needs outlier removal'};
    spikes.params.display.label_colors = [ .7 .7 .7; .3 .8 .3; .3 .3 .8;  .8 .3 .3; .7 .7 .3];
    
    % figure layout
    spikes.params.display.default_figure_size = [.05 .1 .9 .8];  % location of default figure
    spikes.params.display.figure_font_size = 8;  % default font size 
    spikes.params.display.initial_split_figure_panels = 4;
   
    % tool colors
    spikes.params.display.merge_fig_color = [ .7 .8 .7];  % main splitmerge_tool screen 
    spikes.params.display.split_fig_color = [ .8 .7 .7];  % split_tool
    spikes.params.display.outlier_fig_color = [ .7 .7 .8]; % outlier_tool
  
    % axes/panels grid parameters
    spikes.params.display.margin = 140;                 % number of pixels between plots
    spikes.params.display.outer_margin = 80;            % number of pixels around at figure margin
    spikes.params.display.width = 140;                  % number of pixels in plot width
    spikes.params.display.aspect_ratio = 2/3;           % aspect ratio of plots

        
for j = 1:(length(varargin)/2)
    spikes.params = setfield(spikes.params, varargin{(2*j)-1 }, varargin{2*j} );
end

end