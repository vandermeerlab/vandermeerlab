function [k,d] = CORE_mindist(x,y)
%CORE_MINDIST      Core computational routine for MINDIST.
%   [K,D] = CORE_MINDIST(X,Y), given an M x P matrix X and an 
%   N x P matrix Y, returns the M x 1 vectors K and D such that 
%   K(i) = argmin_j ||X(i,:) - Y(j,:)||^2   and
%   D(i) = min_j ||X(i,:) - Y(j,:)||.
%   (||.|| means the Euclidean 2-norm -- note that the values in
%   D are _squared_ distances.)
% 
%   CONDITIONS
%   ----------
%   X and Y must be REAL 2-D arrays of type DOUBLE.
%   X and Y must have the same number of columns.

%%%%%%%%%%%%%%%%%%%%%%%%%%% Prep inputs %%%%%%%%%%%%%%%%%%%%%%%%%%
normsqrX = sum(x.^2,2);
normsqrY = sum(y.^2,2);

[N,P1] = size(x);   % not going to check P1==P2 out of stubbornness --
[M,P2] = size(y);   %   CORE_ functions are described as not doing error checking


%%%%%%%%%%%%%%%%%%%%%% Nearest Vector Search %%%%%%%%%%%%%%%%%%%%%
% Warning: May be verry slow ...
k = zeros(N,1);  d(i) = repmat(Inf, N, 1);
for i = 1:N
	for j = 1:M
		distsq = sum((x(i,:) - y(j,:)).^2);
		if (distsq < d(i)), 
			d(i) = distsq;  k(i) = j;
		end
	end
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1000;
% X = randn(N,10);
% tic;  [K,D] = CORE_mindist(X,X);
% printf('\nCORE_mindist took %5.3f sec.', toc);
% if (K ~= [1:N]' | ~all(abs(D)<2*eps)), printf('  ... but did not find the correct answer.');  end

