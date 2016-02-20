function dx = dxdt(tvec,x,varargin)

% dx = dxdt(tvec,x,varargin)
% window = 1; % seconds
% postSmoothing = 0.5; % seconds --- 0 means don't
%
% Based on Janabi-Sharifi/Hayward/Chen, Discrete-time adaptive windowing
% for velocity estimation, IEEE Transactions on Control Systems Technology,
% (2000) 8(6):1003-1009.
% But modified extensively.  Basic algorithm is to allow windows from 3
% steps to nW = window/DT steps.  For each window, let dx = x(i+nW) - x(i).
% Select the window with the smallest MSE = sum_k=1..nW (x(i+k) - linear-fit(x(i+k) given slope from dx)).^2.
%
% postSmoothing does a convolution of normalized ones(nPS/DT)
%
% ADR 2003
% MvdM edit 2014-07-29: remove use of tsd objects

window = 1; % seconds
postSmoothing = 0.5; % seconds --- 0 means don't
verbose = 1;
extract_varargin;

dT = nanmean(diff(tvec));

nW = ceil(window/dT);
nX = length(x);

if isvectord(x) == 2
   x = x'; 
end

MSE = zeros(nX, nW);
b = zeros(nX,nW);

MSE(:,1:2) = Inf;
nanvector = nan(nW,1);

for iN = 3:nW
	if verbose; fprintf(2,'.'); end
	b(:,iN) = ([nanvector(1:iN); x(1:(end-iN))] - x)/iN;
	for iK = 1:iN
		q = ([nanvector(1:iK); x(1:(end-iK))] - x + b(:,iN) * iK);
		MSE(:,iN) = MSE(:,iN) + q.*q;		
	end
	MSE(:,iN) = MSE(:,iN)/iN;	
end
if verbose; fprintf(2, '!'); end;

[~, nSelect] = min(MSE,[],2);
dx = nan .* ones(size(x));
for iX = 1:nX
	dx(iX) = b(iX,nSelect(iX)) / dT;
end

if postSmoothing ~= 0
	nS = ceil(postSmoothing/dT);
	dx = conv2(dx,ones(nS)/nS,'same');
end

dx = dx'; 

if verbose; fprintf('\n'); end
end