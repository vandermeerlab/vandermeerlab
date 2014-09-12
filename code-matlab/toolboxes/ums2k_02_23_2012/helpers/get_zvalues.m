function [z,r] = get_zvalues( w, covar, r )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% get_zvalues - calculate normalized mahalanobis distance of spikes
%
% Usage:
%    [z,dof] = get_zvalues( w, covar, r )
%
% Description:
%     Calculates the Mahalanobis distance of each waveform in w from the
% mean waveform using the weighting matrix "covar".  Only the top r
% dimensions are used in order to speed up the calculation.  If r is not
% specified then the top 50% of non-singular dimensions will be used. I 
% refer to this function as a z-score since this function calculates a 
% z-score if the waveforms are 1-dimensional and covar is the variance of 
% the waveforms.
%
% Input: 
%   w      - [N x D X C] waveform matrix
%   covar  - covariance matrix to calculate distance
%
% Optional input:
%   r      - rank of used covariance matrix using top principal components 
%
% Output:
%   z - [1 x N] normalized Mahalanolbis distance for each waveform
%   r - value of r used, for when r is not specified 
%

    % if r is not specified, use top 50% of non-singular dimensions
    if nargin < 3
        r = round( rank( covar ) /2 );
    end
   
    % get princiapl components
    num_dims = size(w(:,:),2);
    num_spikes = size(w,1);    
    last = [1:r] + num_dims  - r;    
    [v,d] = eig(covar);            % get PCs
    for j = 1:num_dims, v(:,j) = v(:,j); end
    v = v(:,last);                 % use last r dimensions
    w = detrend(w,'constant');     % mean subtract
    w = (w*v);                     % project on to PCs
 
    
    % get Mahalanobis distance
    z = zeros([1 num_spikes]);
    dinv = inv(d(last,last));
    for j = 1:num_spikes
        z(j) = w(j,:)*dinv*w(j,:)';
    end
            
end