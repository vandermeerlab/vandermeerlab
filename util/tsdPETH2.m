function [P,X0D,SD] = tsdPETH2(X, t, varargin)

% [P,X0D,SD] = tsdPETH(X, t, varargin)
%
% input: 
%        TSD X
%        Event times t
% 
% parameters:
%   window = [-2 5]
%   dt = 0.5; % time steps
%
%  MvdM fix for length(t) == 1, Feb 09

debug = 0;
window = [-2 5];
dt = 0.25;
interp = 'linear';
Extract_varargin;

x = window(1):dt:window(2);

t = t(~isnan(t));

if size(Data(X),2)==1, D=1; else D=2; end;
nT = length(t);

switch D
    case 1
        X0D = nan(nT,length(x));
    case 2
        X0D = nan(nT,length(x),size(Data(X,2))); % probably not correct
end

tGood = 0;

for iT = nT:-1:1
	
	%if (debug && (mod(iT,100) == 0))
	%	disp(sprintf('DEBUG - tsdPETH.m: Spike %d of %d',iT,nT));
	%end
    
    %DisplayProgress(nT-iT, nT);
    X0 = Restrict(X, t(iT)+window(1)-1, t(iT)+window(2)+1);
    X0T = Range(X0,'sec') - StartTime(X0) + window(1)-1;
    if length(X0T) > 1
        tGood = tGood + 1;
        switch D
            case 1
                X0D(iT,:) = interp1(X0T, Data(X0), x, interp);
            case 2
                X0D(iT,:,:) = interp1(X0T, Data(X0), x, interp);
            otherwise
                error('Not yet implemented');
        end
    end
end


P = tsd(x, squeeze(nanmean(X0D))');
SD = tsd(x, squeeze(nanstd(X0D))');


if nargout == 0
	% display
	clf; 
	subplot(2,1,1);
	imagesc(x, 1:nT, X0D);
	axis xy
	xlabel('Peri-event (sec)');
	ylabel('Event #');
	
	subplot(2,1,2);
	m = nanmean(X0D); 
	se =  nanstd(X0D)/sqrt(nT+1);
	plot(x,m,'b',x,m+se,'r:',x,m-se,'r:');
	set(gca, 'XLim', window);
	ylabel('value')
	xlabel('peri-event (sec)');

end