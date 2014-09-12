function spikes = ss_aggregate(spikes, reintegrate_outliers)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_aggregate - Heirarchical cluster aggregation.
%
% Usage:
%      spikes = ss_aggregate(spikes, reintegrate_outliers)
%
% Description:  
% SS_AGGREGATE Heirarchical cluster aggregation.
%     SPIKES = SS_AGGREGATE(SPIKES) takes and returns a spike-sorting object
%     SPIKES after aggregating clusters (requires a previous overclustering
%     and an interface energy matrix calculation).  The aggregation tree is
%     stored in SPIKES.INFO.TREE and the new assignments are stored in
%     SPIKES.ASSIGNS.
%
%     The algorithm computes a similarity/connection matrix using the interface
%     energy.  It then chooses the cluster pair with the highest connection
%     strength, aggregates them, recalculates connection strengths, and then
%     repeats the process.  Aggregation stops when all remaining pairs have
%     a connection strength below a (heuristic) cutoff given by 
%     SPIKES.PARAMS.AGG_CUTOFF.
%
%     The SPIKES.INFO.TREE output is a matrix describing the aggregation.
%     Each aggregation step entry produces a row, listed in the order they
%     were performed.  The first two columns are the indices of the clusters
%     that were aggregated; the index assigned to the aggregate for future
%     entries is the lower of the two original indices.  
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88
%

debug = 1;
starttime = clock;
cutoff = spikes.params.agg_cutoff;   % arbitrarily stop aggregation when overlap density is less than the given % of main cluster density

%%%%% ARGUMENT CHECKING
if (~isfield(spikes.info, 'kmeans') | ~isfield(spikes.info, 'interface_energy'))
    error('SS:energy_not_computed', 'An energy matrix must be computed before aggregation.');
elseif (~isfield(spikes.info, 'kmeans') | ~isfield(spikes.info.kmeans, 'assigns'))
    error('SS:overcluster_not_computed', 'The data must be overclustered before aggregation.');
elseif (isfield(spikes, 'spiketimes') & ~isfield(spikes.params, 'Fs'))
    error('SS:bad_stop_condition', 'A sampling frequency Fs must be supplied as spikes.params.Fs for an ISI stop condition.');
end
if (isfield(spikes, 'assigns') && any(spikes.assigns == 0))
    error('SS:aggregate_after_outliers', 'Aggregation can not be performed after outliers are reintegrated into the data.');
end
if (nargin < 2),  reintegrate_outliers = 1;  end

%%%%% INITIALIZE A FEW THINGS
assignments = double(spikes.info.kmeans.assigns);
interface_energy = spikes.info.interface_energy;
numclusts = max(assignments);
numpts = full(sparse(assignments, 1, 1, numclusts, 1));
tree = [];
untested = 3*ones(numclusts);    % they're all initially untested

handle_fig = figure;
handle_img = imagesc(untested);  axis square;
colormap([0 0 0; 0.9 0.4 0.4; 1 1 0; 0.4 0.8 0.2]);  % [bk, rd, yl, gn] => (0 combined, 1 not allowed, 2 testing, 3 untested)
%%%%% AGGREGATE HIGHEST CONNECTION STRENGTHS UNTIL ALL TRIED
while (any(any(triu(untested,1)==3)))   % only the upper triangular half is meaningful
    % compute connection strengths from interface energies
    %   first, normalize energy:
    normalize = ((numpts * numpts') - diag(numpts));    % Off diag: Na*Nb, On diag: Na^2-Na ...
    normalize = normalize - diag(0.5*diag(normalize));  % ... and divide diagonal by 2
    norm_energy = interface_energy ./ normalize;

    %   then, compute connection strength
    self = repmat(diag(norm_energy), [1,numclusts]);
    connect_strength = 2 .* norm_energy ./ (self + self');
    connect_strength = connect_strength .* (1-eye(numclusts));  % diag entries <- 0, so we won't agg clusters with themselves

    % Find best remaining pair
    remaining = ((untested == 3) .* connect_strength);
    best = max(remaining(:));           % highest untested connection strength
    
    if (best < cutoff)   % No point continuing if connection strengths have gotten really lousy
        break;
    end
    
    [clust1 clust2] = find(connect_strength == best);  % who're the lucky winners?
    first = min(clust1(1),clust2(1));   % if we get 2 best pairs, just take the 1st
    second = max(clust1(1),clust2(1)); 
    untested(first,second) = 2;    untested(first,second) = 2;     % mark that we're trying this one
    set(handle_img, 'CData', triu(untested)); title(['Trying ' num2str(first) ' and ' num2str(second)]); drawnow;
    
     scores = [0 0 0];
  
        % Aggregation subsumes the higher index cluster into the lower.  Start by adding
        % (denormalized) interaction energies for the second (higher index) cluster
        % to those of the first and zeroing the old entries of the second.  Because of the
        % triangular structure of the connection matrix, repeat for both rows and columns ...
        interface_energy(first,:) = interface_energy(first,:) + interface_energy(second,:);
        interface_energy(second,:) = 0;
        interface_energy(:,first) = interface_energy(:,first) + interface_energy(:,second);
        interface_energy(:,second) = 0;
        interface_energy(second,second) = 1;  % keep self-energy at 1 (we may divide by it later)
        % since we added rows & columns, some energy values will have spilled over into the
        % lower half of the energy matrix (which must be upper triangular).  The next 2 steps
        % recover those values.
        overflow = tril(interface_energy, -1);   % spillover below diagonal
        interface_energy = interface_energy + overflow' - overflow;  % reflect above diagonal

        % update counts vector
        numpts(first) = numpts(first) + numpts(second);
        numpts(second) = 2;   % leaving this as 2 prevents div by zero during normalization above
        
        % make a tree entry for the aggregation we just performed
        tree = cat(1, tree, [first, second]);

        % Now actually change the numbers
        assignments(find(assignments == second)) = first;
        
        % Finally, indicate that potential aggregations between the new cluster and 
        % other (nonempty) clusters are untested while pairs involving clusters that
        % have already been emptied should not be tested.
        untested(first,:) = 3;               untested(:,first) = 3;
        untested(tree(:,2),:) = 0;           untested(:,tree(:,2)) = 0;
end
close(handle_fig);

spikes.info.tree = tree;
spikes.assigns = assignments;

% assign colors and fill in labels for all clusters
spikes.info.kmeans.colors = emphasize_cluster_colors(spikes);
spikes.labels = [ sort(unique(spikes.assigns) )' ];
spikes.labels(:,2) = 1; % labels default to 1
