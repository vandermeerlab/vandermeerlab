function dx = dxdt(x,varargin)

% dx = dxdt(x,varargin)
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

window = 1; % seconds
postSmoothing = 0.5; % seconds --- 0 means don't
extract_varargin;

x = ctsd(x);
xD = Data(x);
dT = DT(x);

nW = ceil(window/DT(x));
nX = length(xD);

MSE = zeros(nX, nW);
b = zeros(nX,nW);

MSE(:,1:2) = Inf;
nanvector = repmat(nan, nW,1);

for iN = 3:nW
	fprintf(2,'.');
	b(:,iN) = ([nanvector(1:iN); xD(1:(end-iN))] - xD)/iN;
	for iK = 1:iN
		q = ([nanvector(1:iK); xD(1:(end-iK))] - xD + b(:,iN) * iK);
		MSE(:,iN) = MSE(:,iN) + q.*q;		
	end
	MSE(:,iN) = MSE(:,iN)/iN;	
end
fprintf(2, '!');

[ignore, nSelect] = min(MSE,[],2);
dx = nan .* ones(size(xD));
for iX = 1:nX
	dx(iX) = b(iX,nSelect(iX)) / dT;
end

if postSmoothing
	nS = ceil(postSmoothing/DT(x));
	dx = conv2(dx,ones(nS)/nS,'same');
end
	
dx = tsd(Range(x),dx);
