function spikes = ss_energy(spikes)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_energy - Interface energy based cluster similarity computation.
%
% Usage:
%      spikes = ss_energy(spikes)
%
% Description:  
%     SPIKES = SS_ENERGY(SPIKES) adds an interface energy matrix to a
%     spike-sorting object in SPIKES.INFO.INTERFACE_ENERGY.
% 
%     The energy similarity matrix is calculated by applying an exponential
%     decay to all pairwise euclidean distances between waveforms from two
%     clusters (or within a single cluster for intra-cluster energy) and 
%     summing these distances.
%
%     The calculation ignores the energy due to the zero distance between
%     points and themselves; this removes a dependence of the density on
%     the absolute size of the cluster.  As a result, singleton clusters
%     do not have a well-defined energy and will cause an error.
%  
%     When each entry is normalized by the number of distinct contributing
%     pairs (Na*Nb for off diagonal entries and Na*(Na-1)/2 on the diagonal),
%     it approximates the fraction of pairs in a given cluster whose distance
%     is not much greater than the length constant of the exponential and thus
%     provides an estimate of local density.  This function does not, however,
%     normalize SPIKES.INFO.INTERFACE_ENERGY, since the normalized form is
%     inconvenient during cluster aggregation.  The normalization can readily
%     be done, however, with
%          normalize = ((numpts * numpts') - diag(numpts));
%          normalize = normalize - diag(0.5*diag(normalize));
%          normalized_energy = interface_energy ./ normalize;
%     where 'numpts' is a vector of cluster sizes.
%
%     The unnormalized energy matrix can be updated during aggregation without
%     the need to recompute it from scratch.  The intra-cluster energy E(AB,AB)
%     of a cluster AB formed by aggregating clusters A and B is given by
%              E(AB,AB) = E(A,A) + E(B,B) + E(A,B)
%     and the inter-cluster energy between any cluster C and an aggregate AB is
%                 E(AB,C) = E(A,C) + E(B,C)
%
%      Note that this function also renumbers the miniclusters based on 
%      connection strength.
%
% References:
%     Fee MS et al (1996).  J. Neurosci Methods (69): 175-88
%

starttime = clock;

%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'waveforms') | (size(spikes.waveforms, 1) < 1))
    error('SS:waveforms_undefined', 'The SS object does not contain any waveforms!');
elseif (~isfield(spikes.info, 'kmeans'))
    error('SS:overcluster_not_computed', 'The data must be overclustered before computing energy');
end
numclusts = length(unique(spikes.info.kmeans.assigns));
d = diag(spikes.info.pca.s);
r = find( cumsum(d)/sum(d) >.95,1);
waves = spikes.waveforms(:,:) * spikes.info.pca.v(:,1:r);

%%%%% PREPARE SOME INFORMATION
normsqr = sum(waves.^2,2);
pts = cell(numclusts,1);    % collect spike indices for each cluster
for clust = 1:numclusts   
    pts{clust} = find(spikes.info.kmeans.assigns == clust);
end
numpts = cellfun('length', pts);
if (any(numpts < 2))
    error('SS:energy_ill_defined', 'Clusters with fewer than 2 points do not have a defined energy.');
end

%%%%% HEURISTIC DISTANCE SCALE that seems to work.  The calculation is not too sensitive to this parameter.
scale = sqrt(sum(diag(spikes.info.kmeans.W))) / 10;

%%%%% PREPARE TO LOOP
total = (numclusts^2 + numclusts) / 2;
k = 1;
progress_bar(0, max(floor(total/100),1), 'Computing Interaction Energies . . .')
interface_energy = zeros(numclusts);

%%%%% PAIRWISE DISTANCES LOOP
assigns = spikes.info.kmeans.assigns;
for clust1 = 1:numclusts
    for clust2 = clust1:numclusts   % clust2 starts at clust1 so we get intra- too 
		dists = pairdist(waves(pts{clust1},:), waves(pts{clust2},:));
        interface_energy(clust1,clust2) = sum(exp(-dists(:)/scale));
 
        k = k + 1;
        progress_bar(k/total);
    end
end

%%%%% CORRECTION TERMS
% The energy matrix so far includes a contribution in the intra-cluster
% energies that is not found in the inter-cluster energies; namely, the
% computation of   sum_(all x) sum_(all y) e^(-dist/scale)   for
% intra-cluster energy includes cases where x == y (so dist == 0).
interface_energy = interface_energy - diag(numpts);     % So subtract this out.

% Also, we've double counted pairs in the intra-energy case, since dist(a,b)
% and dist(b,a) are not treated as distinct;
interface_energy = interface_energy - diag(0.5*diag(interface_energy));

%%%%% FINISH UP
spikes.info.interface_energy = interface_energy;

% Now re-number all clusters based on connection strength similarities
assignments = double(spikes.info.kmeans.assigns);
interface_energy = spikes.info.interface_energy;
numclusts = max(assignments);
numpts = full(sparse(assignments, 1, 1, numclusts, 1));
normalize = ((numpts * numpts') - diag(numpts));    % Off diag: Na*Nb, On diag: Na^2-Na ...
normalize = normalize - diag(0.5*diag(normalize));  % ... and divide diagonal by 2
norm_energy = interface_energy ./ normalize;
self = repmat(diag(norm_energy), [1,numclusts]);
connect_strength = 2 .* norm_energy ./ (self + self');
connect_strength = connect_strength .* (1-eye(numclusts));  % diag entries <- 0, so we won't agg clusters with themselves

% initialize
[val,pos(1)] = max(max(connect_strength) ) ;
for j = 2:numclusts
    [val1,pos1] = max(connect_strength(:,pos(j-1)) );
    [val2,pos2] = max(connect_strength(pos(j-1),: ));
    if val1 >= val2, pos(j) = pos1; else, pos(j) = pos2; end
    connect_strength(:,pos(j-1)) = -inf;
    connect_strength(pos(j-1),:) = -inf;
end

% update kmeans assignments
assigns = spikes.info.kmeans.assigns;
ie = spikes.info.interface_energy;
c = spikes.info.kmeans.centroids;
for j = 1:numclusts
   assigns( spikes.info.kmeans.assigns == pos(j) ) = j;
   for k = 1:numclusts
       if k > j
        if pos(k) > pos(j)
          ie( j,k) = spikes.info.interface_energy( pos(j),pos(k) );
        else
          ie( j,k) = spikes.info.interface_energy( pos(k),pos(j) );
        end
       end
   end
   c(j,:) = spikes.info.kmeans.centroids(pos(j),:);
end
spikes.info.kmeans.assigns = assigns;
spikes.info.interface_energy = ie;
spikes.info.kmeans.centroids = c;

