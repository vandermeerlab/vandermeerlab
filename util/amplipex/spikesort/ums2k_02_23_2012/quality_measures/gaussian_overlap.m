function confusion = gaussian_overlap( w1, w2 )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% gaussian_overlap - estimate cluster overlap from a 2-mean Gaussian mixture model
%
% Usage:
%     C = gaussian_overlap( waveforms1, waveforms2 )
%
% Description:  
%   Estimates the overlap between 2 spike clusters by fitting with two
% multivariate Gaussians.  Implementation makes use of MATLAB statistics
% toolbox function "gmdistribution.fit". 
%
% The percent of false positive and false negative errors are estimated for 
% both classes and stored as a confusion matrix. Error rates are calculated 
% by integrating the posterior probability of a misclassification.  The 
% integral is then normalized by the number of events in the cluster of
% interest. See description of confusion matrix below.
%
% NOTE: The dimensionality of the data set is reduced to the top 99% of 
% principal components to increase the time efficiency of the fitting
% algorithm.
%
% Input: 
%   waveforms1  - [Event x Sample ] waveforms of 1st cluster
%   waveforms2  - [Event x Sample ] waveforms of 2nd cluster
%      If your waveform data is in the form [Event X Sample X Channels ], 
%      you should call this function as
%            confusion = gaussian_overlap( w1(:,:), w2(:,:) );
%      This will concatenate the waveforms from different channels.
%      
%
% Output:
%   C   - a confusion matrix
%   C(1,1) - False positive fraction in cluster 1 (waveforms of neuron 2 that were assigned to neuron 1)
%   C(1,2) - False negative fraction in cluster 1 (waveforms of neuron 1 that were assigned to neuron 2)
%   C(2,1) - False negative fraction in cluster 2 
%   C(2,2) - False positive fraction in cluster 2
%

    % reduce dimensionality to 98% of top Principal Components
    cutoff = 98;
    N1 = size(w1,1);
    N2 = size(w2,1);
    disp('Reducing data dimensionality ...');
   [proj,u,s,v] = pcasvd([w1; w2]);
   for j = 1:size(s,2)
       vals(j) = s(j,j);
   end
   cumvals = cumsum(vals);
   thresh  = sum(vals) * (cutoff/100);
   num_dims = find( cumvals < thresh, 1, 'last' ) - 1;
   if isempty(num_dims), num_dims = size( w1(:,:), 2 ); end
   w1 = proj(1:N1,1:num_dims);
   w2 = proj([N1+1:end],1:num_dims);
  
    % fit 2 multivariate gaussians, use observed parameters to initialize
    params.mu = [ mean(w1(:,:),1); mean(w2(:,:),1)];
    params.Sigma(:,:,1) = cov(w1(:,:));
    params.Sigma(:,:,2) = cov(w2(:,:));
    params.PComponets = [N1 N2] / (N1 + N2);
    disp('Fitting 2-Gaussian Mixture Model ...')
    gmfit = gmdistribution.fit([w1(:,:); w2(:,:)],2,'Start',params); 
   
    % get posteriors
    disp('Calculating Confusion matrix ...')
    pr1 = gmfit.posterior(w1);
    pr2 = gmfit.posterior(w2);

    % in the unlikely case that the cluster identities were flipped during the fitting procedure, flip them back
    if mean(pr1(:,1)) + mean(pr2(:,2)) < 1
        pr1 = pr1(:,[2 1]);
        pr2 = pr2(:,[2 1]);
    end
    
    % create confusion matrix
    confusion(1,1) = mean(pr1(:,2));   % probability that a member of 1 is false
    confusion(1,2) = sum(pr2(:,1))/N1 ; % relative proportion of spikes that were placed in cluster 2 by mistake
    confusion(2,2) = mean(pr2(:,1)); % probability that a member of 2 was really from 1
    confusion(2,1) = sum(pr1(:,2))/N2; % relative proportion of spikes that were placed in cluster 1 by mistake

    % 
    disp('Finished.')

end

