function [dF,baseline] = baselineFP(FP,interpType,fitType,basePrc,winSize,winOv,Fs)
%baselinePhotometry - Baseline adjust photometry signal to get dF/F
%
%   [dF_F,varargout] = baselineFP(FP,interpType,fitType,basePrc,winSize,winOv,Fs)
%
%   Description: This code will baseline adjust the photometry signal using
%   a moving window and finding values within a specified percentile
%
%   Input:
%   - FP - Photometry signal to baseline
%   - interType - Interpolation method to use: 'linear' 'spline' etc
%   - fitType - Options: 'interp' --> interpolated line
%       - 'exp' --> Models interpolated line as exponentials
%       - 'line' --> models interpolated line as first degree polynomial
%   - basePrc - Percentile value to use when finding baseline points for
%   interpolation
%   - winSize - Window size in seconds for finding baseline
%   - winOv - Overlap size in seconds for finding baseline
%   - Fs - Sampling Rate
%
%   Output:
%   - dF_F - Baseline adjusted trace (%dF_F)
%   - varargout - Optional fitline output
%
%   Author: Pratik Mistry, 2019
%

%Ensure the FP vector is a column vector instead of row vector --> Faster
%computation
if size(FP,1) == 1
    FP = FP';
end

winSize = winSize * Fs; %Convert window size from seconds to samples
winOv = winOv * Fs; %Convert overlap window from seconds to samples
Ls = length(FP); %Get length of photometry trace
L = 1:Ls; L = L'; %Create a column vector from 1 to total data points in trace
nPts = floor(Ls/(winSize-winOv)); %Determine number of baseline points to be found

%X is a vector of positional points in the data that we will be gathering baseline points
%Y is a vector of zeros that will contain calculated baseline points
X = L(1:ceil(Ls/nPts):end); Y = zeros(nPts,1);

%Determine the step size of the window:
%If the overlap is 0 or empty then it will use the window size as the step
%size. If the overlap is greater than 0 the step size will be the window
%size subtracted by the overlap size
if winOv == 0 || isempty(winOv)
    winStep = winSize;
else
    winStep = winSize - winOv;
end

%The following for loop goes through the photometry vector and finds
%baseline values of the windowed photometry trace according to a certain
%percentile
for n = 0:nPts-1
    I1 = (n*winStep)+1;
    I2 = I1 + winSize;
    if I2>Ls
        I2 = Ls;
    end
    Y(n+1) = prctile(FP(I1:I2),basePrc);
end

interpFit = interp1(X,Y,L,interpType,'extrap'); %Create an interpolated line using the previously calculated baseline points


%The following switch statement will adjust the baseline according to user
%input
switch fitType
    case 'interp'
        baseline = interpFit; %Use interpolated line
    case 'exp'
        expFit = fit(L,interpFit,'exp2'); %Baseline is now an exponential line fit to interpolated line
        baseline = double(expFit(L));
    case 'line'
        lineFit = fit(L,interpFit,'poly1'); %Baseline is a linear fit of the interpolated line
        baseline = double(lineFit(L));
    otherwise
        
end

dF = (FP-baseline)./baseline;
dF = dF*100;

end
