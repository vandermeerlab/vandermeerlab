function [outputS, outputT, outputG] = spikePETH(S,t,varargin)

% spikePETH(S, t, varargin)
%
% input: 
%        TS S
%        Event times t
% 
% if no outputs, plots
%
% parameters:
%   window = [-2 5]
%   dt = 0.001;

window = [-2 5];
dt = 0.00025;
excessBounds = 1;
outputGrid = 0;
extract_varargin;

nT = length(t);

outputS = [];
outputT = [];
%outputID = repmat(inf, nT, diff(window)/dt+1);
outputIT = linspace(window(1), window(2), diff(window)/dt+1);

if outputGrid
    xbin = window(1):dt:window(2);
    outputG = zeros(nT,length(xbin)-1);
end

for iT = 1:nT
	
	%DisplayProgress(iT, nT);
	S0 = Restrict(S, t(iT)+window(1)-excessBounds, t(iT)+window(2)+excessBounds);
	if length(Data(S0)) > 0
		%I0 = ISIctsd(S0,'dt',dt,'t0',t(iT)+window(1),'t1',t(iT)+window(2));
		S0 = Restrict(S0, t(iT)+window(1), t(iT)+window(2));
		
		%outputID(iT,:) = Data(I0)';
		
		outputT = [outputT; repmat(iT, length(Data(S0)),1)];
		outputS = [outputS; Data(S0)-t(iT)];
        
        if outputGrid
           
            temp = histc(Data(S0)-t(iT),xbin); temp = temp(1:end-1);
            if ~isempty(temp)
                outputG(iT,:) = temp;
            end
            
        end
        
	end
end

if nargout == 0
	% display
	clf
	
	subplot(2,1,1);
% 	imagesc(window,[1 nT], outputID); 
% 	colormap(1-0.25*gray);
% 	hold on;
	plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
	xlabel('peri-event (sec)');
	ylabel('Event #');
	
	%
	subplot(2,1,2);
	m = histc(outputS, outputIT);
	bar(outputIT,m/dt/length(t));
% 	x = outputIT;
% 	m = nanmean(1./outputID); 
% 	se =  nanstd(1./outputID)/sqrt(nT+1);
% 	plot(x,m,'b',x,m+se,'r:',x,m-se,'r:');
	set(gca, 'XLim', window);
	ylabel('FR (Hz)')
	xlabel('peri-event (sec)');
end