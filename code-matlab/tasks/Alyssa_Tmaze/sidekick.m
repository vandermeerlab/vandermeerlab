function sidekick(X,varargin)
%SIDEKICK Add plots to existing figure (particularly MultiRaster plotmode 7).
%   Add a 2-D line plot of the data in each arg versus the corresponding
%   values in X. X must be the same length as the one used in the existing 
%   figure. Can also plot additional iv inputs.
%
%   Varargins can be variables or expressions, but unless they're variables
%   or have a name field, the legend won't work. 
%
%   If X and arg are both vectors, they must have the same length.
%
%   If arg is an iv data type, then sidekick plots it as vertical bars. Add
%   a name field for the legend (evt.name = 'getSWR').
%
%   Ex syntax: sidekick(csc.tvec,Zscore,swrscore,evt)
%   Ex legend fail: sidekick(csc.tvec,Zscore,sin(csc.tvec))
%
% Colour order is blue, red, green, orange, purple, turqoise, gold, majenta
% cyan, sky blue, cobalt green, brown.  
%
%
% RESCALING NOT IMPLEMENTED YET
% ACarey, Jan 2015

%%

colVal = {[0 0 1],[1 0 0],[0 1 0],[1 165/255 0],[191/255 62/255 1],[0 206/255 209/255],[1 215/255 0],[1 0 1],[0 1 1],[61/255 145/255 64/255],[139/255 69/255 19/255]}; % these are the RGB values
lgndcol = []; 
%fig = gcf;
hold on;
for iArg = 1:length(varargin)
    Y = varargin{iArg};
    % check if Y is an iv data type
    if isfield(Y,'tstart')
       ylims = get(gca,'Ylim');
       plot([Y.tstart Y.tstart],ylims,'Color',colVal{iArg})
       plot([Y.tend Y.tend],ylims,'Color',colVal{iArg})
       if isfield(Y,'name')
           lgnd{iArg} = Y.name;
           lgndcol = [lgndcol; [colVal{iArg}]];
       end
    else
        plot(X,Y,'Color',colVal{iArg})
        lgnd{iArg} = inputname(iArg+1);
        lgndcol = [lgndcol; [colVal{iArg}]];
    end
end
% now do legend
lgndcol = flipud(lgndcol); % because MATLAB legend order is backwards to the plot order, the colors are wrong without this
leg = legend(lgnd);
legtxt = findobj(leg,'type','text');
%assignin('base','legtxt',legtxt)
for itxt = 1:length(legtxt) 
    set(legtxt(itxt),'color',lgndcol(itxt,:));
end

hold off;
