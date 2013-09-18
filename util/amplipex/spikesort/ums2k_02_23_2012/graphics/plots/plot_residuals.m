function plot_residuals( spikes, show )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_residuals - plots standard deviation of waveforms for a group of events
%
% Usage:
%      plot_residuals( spikes, show )
%
% Description:  
%   For the specified spikes, plot the standard deviation of the waveform
% as a function of sample.  Red vertical bars separate samples from 
% different waveforms. Black dotted horizontal line indicates background
% noise level on that channel.  Pink region is 95% confidence interval
% for significant difference from background noise.
%
% Inputs:
%   spikes - a spikes object, must have cluster assignments
%
% Optional input:
%   show          - array describing which events to show in plot
%                 - see get_spike_indices.m, (default = 'all')
%

    % handle arguments
    if nargin < 2, show = 'all'; end

    % calcualte standard deviation
    clus = get_spike_indices(spikes, show );
    memberwaves = spikes.waveforms(clus,:); 
    num_channels = size(spikes.waveforms,3);
    num_samples  = size(spikes.waveforms,2);
    s = std( memberwaves );  

    % plot it
    cla reset
    set( line(1:length(s),s) ,'Color', [.3 .5 .3],'LineWidth',1.5)
   
   % plot predicted standard deviation from background noise along with 95% confidence interval
    for j = 1:num_channels
        x = [1 num_samples+1] + (j-1)*num_samples;
        stdev = abs( spikes.info.detect.stds(j) );
        [lb,ub] = std_bounds(stdev, length(clus), .95);
       l(j) =  line( x, [1 1] * stdev );
       
       p(j)  = patch( [ x fliplr(x)], [lb lb ub ub],[1 .8 .8]);
    end
  
    axis tight;          
    set(gca,'YLim',[0  max( max(get(gca,'YLim'),2*max(spikes.info.detect.stds)))]);
    
    % draw dividers
    for j = 1:num_channels-1
        set( line( 1 + num_samples * j * [1 1], get(gca,'YLim')), 'Color',[1 0 0] ) % electrode dividers
    end
    
    set( l, 'Color',[ 0 0 0],'LineStyle',':','LineWidth',2)
    xlabel('Sample')
    ylabel('Residuals')
    uistack(p,'bottom')
  
% calculate significance cutoff for standard deviation estimates via the chi2 distribution
function [lb,ub] = std_bounds( stdev, N, p )

    if nargin < 3, p = .95; end

    bounds = (1 + p*[-1 1]) / 2;
     
    x = chi2inv( bounds,N-1);
    x = (stdev^2) * x / (N-1);

    lb = x(1)^.5; ub = x(2)^.5;
    
    