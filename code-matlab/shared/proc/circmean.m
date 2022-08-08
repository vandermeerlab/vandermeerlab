function [mu,r,ci95,rp] = circmean(phi,varargin)
% [mu,r,ci95,rp] = circmean(phi)
%
% INPUTS
%    phi - nSamp x n array of samples (defaults to rad)
%
% OUTPUTS
%    mu - 1 x n array circular mean of samples in phi
%    r - 1 x n array length of vector mean
%
% MvdM 10 added unit modifier; r, ci outputs


units = 'rad'; % {'rad','deg'}
extract_varargin;

if strcmp(units,'rad')
   phi = phi.*(360./(2*pi)); 
end


X = cos(phi * pi/180);
Y = sin(phi * pi/180);
mX = nanmean(X);
mY = nanmean(Y);
mu = atan2(mY,mX) * 180/pi;


r = sqrt(mX.^2+mY.^2);

% significance (from Batchelet 1981)
n = length(phi);
z = n*(r.^2);
rp = exp(-z).*(1+(((2*z)-(z.^2))/(4*n))-(((24.*z)-(132.*z.^2)+(76.*z.^3)-(9.*z.^4))/(288*n^2)));


if nargout == 3 % do ci

    nD = size(phi,2);
    nS = size(phi,1);
    ci95 = nan(size(mu));
    for iD = 1:nD

        m1 = 0; m2 = 0;
        for iS = 1:nS

            m1 = m1 + cosd(diffang(phi(iS,iD),mu(iD)));
            m2 = m2 + cosd(2*diffang(phi(iS,iD),mu(iD)));

        end
        m1 = m1/nS; m2 = m2/nS;

        scd = (1-m2)/(2*m1^2);
        cse = sqrt(scd/nS);

        % confidence intervals
        if ((icdf('norm',0.975,0,1)*cse) <= 1)
            ci95(iD) = asind(icdf('norm',0.975,0,1)*cse);
        end


    end


end

if strcmp(units,'rad')
    mu = mu./(360./(2*pi));
    if nargout == 3
        ci95 = ci95./(360./(2*pi));
    end
end

function DAng=diffang(Angle1,Angle2)
% diffang: just subtracts to angle2 from angle 1 and reifies it between 
%				-180 and 180
% DAng=diffang(Angle1,Angle2)
DAng=Angle1-Angle2;
x=find(DAng<-180);
DAng(x)=DAng(x)+360;
x=find(DAng>180);
DAng(x)=DAng(x)-360;