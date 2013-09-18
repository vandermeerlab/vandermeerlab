function p = plot_detection_criterion(spikes,clus)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% plot_detection_criterion - histogram of detection criterion applied to 
%                            each waveform
%
% Usage:
%     [thresh_val,m,stdev,missing] = plot_detection_criterion(spikes,clus)
%
% Description:  
%    Plots a histogram of the detection metric applied to each waveform. A
% normalization is performed so that the threshold has a value of +/- 1 
% (vertical black dotted line). The histogram is overlaid with a fitted
% Gaussian (red line) that is used to estimate the number of undetected
% spikes. This percentage is shown in the axes title.
%
%  See ss_undetected.m and undetected.m for more information.
%
% Input:
%   spikes     - a spike structure
%   clus       - cluster ID or array describing which events to show in plot
%              - see get_spike_indices.m
%
% Output:
%  p            - estimate of probability that a spike is missing because it didn't reach threshhold
%

    % check arguments
    if ~isfield(spikes,'waveforms'), error('No waveforms found in spikes object.'); end 
    if nargin < 3, display = 1; end
    if nargin < 2, show = 'all'; end 
  
    % grab data
    [p,mu,stdev,n,x] = ss_undetected(spikes,clus);
    
    % determine global extreme if there are other detection criterion plots on the current figure
    my_sign = sign(mu);
    ax = findobj( gcf, 'Tag', 'detection_criterion' );
    my_ex = max(abs(x));
    if ~isempty(ax)
        xlims = get(ax(1),'XLim');
        if max(abs(xlims)) > my_ex
            global_ex = max(abs(xlims))*my_sign;
        else
            global_ex = my_ex*my_sign;
            set(ax, 'XLim', sort( [global_ex 0]) )
        end
    else
        global_ex = my_ex*my_sign;
    end

    % Now make an estimate of how many spikes are missing, given the Gaussian and the cutoff
    N = sum(n) / (1-p);
    if my_sign == -1
      a = linspace(global_ex,0,200);
    else
      a = linspace(0,global_ex,200);
    end
    b = normpdf(a,mu,stdev);
    b = (b/sum(b))*N*abs((x(2)-x(1))/(a(2)-a(1)));
   
    % plot everything    
    cla reset

    % histogram
    hh = bar(x,n,1.0);
    set(hh,'EdgeColor',[0 0 0 ])
    set( gca,'XLim',sort([ global_ex 0]));
    
    % gaussian fit
    l =line(a,b);
    set(l,'Color',[1 0 0],'LineWidth',1.5)
   
    % threshold line
    l = line([1 1]*my_sign, get(gca,'YLim' ) );
    set(l,'LineStyle','--','Color',[0 0 0],'LineWidth',2)

    % prettify axes
    axis tight
    set( gca,'XLim',sort( [global_ex 0]) );
    set(gca,'Tag','detection_criterion')
    title( ['Estimated missing spikes: ' num2str(p*100,'%2.1f') '%']);
    xlabel('Detection metric')
    ylabel('No. of spikes')

end

