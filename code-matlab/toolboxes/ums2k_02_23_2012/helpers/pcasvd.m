function [proj,u,s,v] = pcasvd(data)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/09/2010
%
% pcasvd - Principal Components Analysis via (mean-subtracted) SVD.
%
% Usage:
%   [proj,u,s,v] = pcasvd(data)
%
% Description:
%   PROJ = PCASVD(DATA), where DATA is an M x N matrix, returns the M x N
%   matrix PROJ where the (M,N)th entry is the projection of the Mth row
%   of DATA onto the Nth eigenvector of the covariance matrix formed from
%   the rows of DATA.
%
%   [PROJ,U,S,V] = PCASVD(DATA) also returns matrices U, S, V such that
%   DATA = U * S * V' and PROJ = DATA * V.
%
%   All of these computations are generally performed taking the mean over
%   all rows of the matrix DATA to be the zero vector.  This is therefore
%   enforced if it is not already the case.


data = detrend(data, 'constant');   % remove mean row
[u,s,v] = svd(data, 0);             % SVD the data matrix
proj = data * v;                    % compute (mean-subtracted) pca projections