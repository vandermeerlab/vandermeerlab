function CSCf = fftFilterLFP(cfg_in,CSC,ncfs)
%FFTFILTERLFP Filter the local field potential using a custom frequency spectrum. 
%   Detailed explanation goes here indeed
%
% ncfs - noise-corrected frequency spectrum (see ducktrap and SWRfreak).
%        These are the frequencies you want to keep.
% send output to photonic
%
% Elyot Grant and aacarey, Nov 2015 (initial version)

%% 
mfun = mfilename;

% initialize defaults, process config
cfg_def.verbose = 1;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

if cfg.verbose % talk to me
    disp([mfun,': filtering the signal'])
end

% rescale and interpolate the ncfs
lenCSC = length(CSC.data); lenNCFS = length(ncfs);
%120 ms window for ducktrap/SWRfreak thing...
scaleFactor = lenCSC/lenNCFS*2;
xi = (1:(lenCSC/2))/scaleFactor;
yi = interp1(1:lenNCFS,ncfs,xi,'spline','extrap');
yi = max(yi,0); % throw out the stuff below 0

% make a filter 
fftfilter = [0 yi yi(end-1:-1:1)];

% do the thing
CSCf = CSC;
CSCf.data = ifft(fftfilter.*fft(CSC.data));

end

