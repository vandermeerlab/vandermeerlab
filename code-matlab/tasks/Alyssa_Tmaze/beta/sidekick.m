function h = sidekick(cfg_in,varargin)
%SIDEKICK Add iv and tsd to existing figure.
%   Add a 2-D line plot of tsd or iv data. tsd plots as tvec vs data, while
%   iv plots as vertical lines at tstart and tend. The order of the
%   varargins does not matter. A corresponding legend is also included, but
%   input names are displayed only if the inputs are variables; if the
%   inputs are expressions, the legend displays 'Input 1' and so on.
%
%   h = SIDEKICK(cfg,varargin);
%
%   h = SIDEKICK(cfg,tsd1,tsd2,iv1,tsd3);
%   h = SIDEKICK(cfg,iv1,tsd1,iv2,iv3,iv4,tsd2,tsd3);
%
%   INPUTS
%       tsd (any number)
%        iv (any number)
%      Note: more data takes longer to plot, and also slows down navigation.
%
%   OUTPUTS
%       h - plot and legend handles
%       h(1) corresponds to input 1 and so on
%       h(end) corresponds to the legend
%
%   CFG OPTIONS
%       cfg.LineWidth = 1.5; 
%       cfg.Location = 'NorthEast'; Location of legend.
%       cfg.ColVal   = [nVararg x 3] array of [r g b] values
%   
%
% RESCALING NOT IMPLEMENTED YET
% aacarey, Jan 2015, Nov 2015

%%

cfg_def.LineWidth = 1.5;
cfg_def.ColVal = linspecer(length(varargin)); % get a bunch of colors to use
cfg_def.Location = 'NorthEast';
cfg_def.openNewFig = 0;
cfg = ProcessConfig2(cfg_def,cfg_in);

if cfg.openNewFig
    figure;
end
hold on;
lgnd = cell(size(varargin));
h = zeros(1,length(varargin)+1);

for iArg = 1:length(varargin)
    temp = varargin{iArg};
    name = inputname(iArg+1); if isempty(name); name = ['Input ',num2str(iArg)]; end; lgnd{iArg} = name;
    
    % check if Y is an iv data type
    if isfield(temp,'tstart')
       ylims = get(gca,'Ylim');
       temp = [temp.tstart; temp.tend];
       x = [temp temp nan(size(temp))]; 
       x = reshape(x',length(temp)*3,1);
       
       y = repmat([ylims nan],length(temp),1);
       y = reshape(y',length(temp)*3,1);
       h(iArg) = plot(x,y,'Color',cfg.ColVal(iArg,:),'LineWidth',cfg.LineWidth);
    else
       h(iArg) = plot(temp.tvec,temp.data,'Color',cfg.ColVal(iArg,:),'LineWidth',cfg.LineWidth);
        
    end
end

% now do legend
h(iArg+1) = legend(h(1:iArg),lgnd);

set(h(iArg+1),'Location',cfg.Location)

hold off;
