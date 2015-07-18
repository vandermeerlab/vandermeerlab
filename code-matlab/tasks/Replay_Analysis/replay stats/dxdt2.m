function dx = dxdt2(cfg_in,tvec,pvec)

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
% youkitan edit 2015-07-13: cfg

cfg_def.window = 1; % seconds
cfg_def.postSmoothing = 0.5; % seconds --- 0 means don't
cfg = ProcessConfig2(cfg_def,cfg_in);

dT = nanmean(diff(tvec));

nW = ceil(cfg.window/dT);
nX = length(pvec);

if isvectord(pvec) == 2
   pvec = pvec'; 
end

MSE = zeros(nX, nW);
b = zeros(nX,nW);

MSE(:,1:2) = Inf;
nanvector = nan(nW,1);

for iN = 3:nW
	fprintf(2,'.');
	b(:,iN) = ([nanvector(1:iN); pvec(1:(end-iN))] - pvec)/iN;
	for iK = 1:iN
		q = ([nanvector(1:iK); pvec(1:(end-iK))] - pvec + b(:,iN) * iK);
		MSE(:,iN) = MSE(:,iN) + q.*q;		
	end
	MSE(:,iN) = MSE(:,iN)/iN;	
end
fprintf(2, '!');

[~, nSelect] = min(MSE,[],2);
dx = nan .* ones(size(pvec));
for iX = 1:nX
	dx(iX) = b(iX,nSelect(iX)) / dT;
end

if cfg.postSmoothing
	nS = ceil(cfg.postSmoothing/dT);
	dx = conv2(dx,ones(nS)/nS,'same');
end

dx = dx'; 

fprintf('\n');

end