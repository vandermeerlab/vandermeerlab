function compare_clusters( spikes, show, alt_assigns )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% compare_clusters -- creates a figure to compare a set of clusters pairwise
%
% Usage:
%       compare_clusters( spikes, clusters, alt_assigns )
%
% Description:  
%    Builds a figure with a grid of plots to compare a given set of clusters. 
%  The plot in the ith row and jth column compares the ith cluster to the
%  jth cluster.  Plots are generated as follows:
%
%     upper-right triangle  =>  plot_fld
%     lower-left triangle   =>  plot_xcorr
%     diagonal              =>  plot_isi
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
if (nargin < 2),     show = sort( unique(spikes.assigns) );  show = setdiff(show,0); end

% set up figure
h = findobj(0,'Name','Separation Analysis');
if isempty(h), h = figure('Units','Normalized','Position',spikes.params.display.default_figure_size); end
clf(h)
drawnow
set(h,'Pointer','watch'),pause(.01)
set(h,'defaultaxesfontsize',spikes.params.display.figure_font_size);
sliderFigure(h,spikes.params.display.outer_margin)

% make a grid of plots --> the plot in row i, column j compares the ith cluster to the jth cluster
show = sort(show);
for row = 1:length(show)
    for col = 1:length(show)

       % make axes at proper position
       pos = get_pos_on_grid( [row col 1 1], spikes.params.display, h, 1 );
       ax(row,col) = axes('Visible','off','Units','pixel','Position',pos);

       
       % in upper right triangle, plot FLD
       if row < col
        plot_fld( spikes, show(row), show(col) );
        legend off;
           set(gca,'Color', [1 .9 .9] )
       
       % in lower left triangle, plot cross-correlation
       elseif row > col
        plot_xcorr( spikes, show(row), show(col) );          
           set(gca,'Color', [.9 .9 1] )
           
       % on diagonal plot ISI distribution
       elseif row == col
           plot_isi( spikes, show(row) )
           
       end

       % label cluster plots along row 1 and column 1
       if row == 1, 
           current_title = get(get(gca,'Title'),'String');
           title( {clust_string( spikes, show(col) ), current_title} ); 
       end

       if col == 1 
           ystr = get(  get(ax(row,col),'YLabel'), 'String');
           if ~iscell(ystr), ystr = {ystr}; end
           ystr = [ clust_string( spikes, show(row) );  ystr];
           set( get(ax(row,col),'YLabel'), 'String', ystr )
       end
       

    end
end
set(ax,'Visible','on')
set(gcf,'Name','Separation Analysis','NumberTitle','off')
set(h,'Pointer','arrow'),pause(.01)

figure(gcf)

% helper function to display cluster name
function str = clust_string( spikes, clus )

num_spikes = sum( spikes.assigns == clus );
if clus == 0
    str = 'Outliers';
else
    str = ['Cluster #' num2str(clus )];
end

str = [str '  (N = ' num2str(num_spikes) ')'];


