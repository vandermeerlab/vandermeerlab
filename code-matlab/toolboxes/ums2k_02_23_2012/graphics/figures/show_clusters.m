function show_clusters(spikes, clusters, alt_assigns)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% show_clusters -- creates a figure to plot various statistics of a set of clusters
%
% Usage:
%       show_clusters( spikes, clusters, alt_assigns )
%
% Description:  
%    Builds a figure with a row of plots for each specified cluster. Plots 
%  are generated as follows:
%     
%    column 1          =>          plot_waveforms
%    column 2          =>          plot_residuals
%    column 3          =>          plot_detection_criterion
%    column 4          =>          plot_isi
%    column 5          =>          plot_stability
%
%  SliderFigure is also called to allow the user to scale plots.
%
% Input: 
%  spikes      - a spikes structure
%
% Optional input:
%  clusters        - list of cluster IDs (default is to use all clusters)
%  alt_assigns - alternate list of spike assignments (default is to use spikes.assigns)
%

% Check arguments

if nargin == 3, spikes.assigns = alt_assigns; end
if ~isfield(spikes,'assigns'), error('No assignments found in spikes object.'); end
if (nargin < 2),     clusters = sort( unique(spikes.assigns) );  end

% set up figure
h = findobj(0,'Name','Show Clusters');
if isempty(h), h = figure('Units','Normalized','Position',spikes.params.display.default_figure_size); end
clf(h)
drawnow
set(h,'Pointer','watch'),pause(.01)
set(h,'defaultaxesfontsize',spikes.params.display.figure_font_size);

% plot row for each cluster
show = sort(clusters);
for row = 1:length(show)

    clus = show(row);

    pos = get_pos_on_grid( [row 1 1 1], spikes.params.display, h, 1 );
    ax(row,1) = axes('Visible','off','Units','pixel','Position',pos);
    plot_waveforms(spikes,clus);
    
    pos = get_pos_on_grid( [row 2 1 1], spikes.params.display, h, 1 );
    ax(row,2) = axes('Visible','off','Units','pixel','Position',pos);
    plot_residuals(spikes,clus);

    pos = get_pos_on_grid( [row 3 1 1], spikes.params.display, h, 1 );
    ax(row,3) = axes('Visible','off','Units','pixel','Position',pos);
    plot_detection_criterion(spikes,clus);

    pos = get_pos_on_grid( [row 4 1 1], spikes.params.display, h, 1 );
    ax(row,4) = axes('Visible','off','Units','pixel','Position',pos);
    plot_isi(spikes,clus);
    
    pos = get_pos_on_grid( [row 5 1 1], spikes.params.display, h, 1 );
    ax(row,5) = axes('Visible','off','Units','pixel','Position',pos);
    [ax(row,5) ax(row,6)] = plot_stability(spikes,clus);
   
     
end

% finish up figure
set(ax,'Visible','on')
sliderFigure(h,spikes.params.display.outer_margin)
set(h,'Name','Show Clusters','NumberTitle','off','Pointer','arrow')
figure(h)

    