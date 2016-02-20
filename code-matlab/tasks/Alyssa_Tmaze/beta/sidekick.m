function h = sidekick(cfg_in,varargin)
%SIDEKICK Add iv and tsd to existing figure.
%   Add a 2-D line plot of tsd or iv data. tsd plots as tvec vs data, while
%   iv plots as vertical lines and patch objects at tstart and tend. The 
%   order of the varargins does not matter. A corresponding legend is also
%   included, but input names are displayed only if the inputs are variables;
%   if the inputs are expressions, the legend displays 'Input 1' and so on.
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
%       cfg.LineWidth  = 1.5; 
%       cfg.showLegend = 1; If 1 makes legend, if 0 doesn't. Useful if
%                        sidekick is called inside of another function
%                        and you don't want a legend and don't want to 
%                        deal with deleting it.
%       cfg.Location   = 'NorthEast'; Location of legend.
%       cfg.Color      = 'linspecer'; What colors to use, a number of 
%                         different input formats:
%                      - string specifying a colormap to use, ex: 'summer'
%                      - [nVararg x 3] array of [r g b] values
%                      - [1 x 3] array of [r g b] values
%                      - string specifying a color to use, ex: 'green' or
%                        'g'.
%       cfg.patch      = 1; If 1, plots IV data as semi-transparent patch
%                        objects; if 0, doesn't.
%
% RESCALING NOT IMPLEMENTED YET
% aacarey, Jan 2015, Nov 2015, Feb 2016

%%

cfg_def.LineWidth = 1.5;
cfg_def.Color = 'linspecer';
cfg_def.Location = 'NorthEast';
cfg_def.showLegend = 1;
cfg_def.patch = 1;

mfun = mfilename;

cfg = ProcessConfig(cfg_def,cfg_in,mfun);

ColorSpecStrings = {'y','m','c','r','g','b','w','k','yellow','magenta','cyan','red','green','blue','white','black'};

% make color matrix
if (ischar(cfg.Color) && ismember(cfg.Color,ColorSpecStrings)) || (isnumeric(cfg.Color) && length(cfg.Color) == 3)
    cfg.ColVal = repmat(cfg.Color,length(varargin),1);
elseif length(cfg.Color) > 1 && ~isnumeric(cfg.Color)
    cfg.ColVal = colormap(eval([cfg.Color '(' num2str(length(varargin)) ')'])); % turn the string into a matrix of RGB values
elseif isnumeric(cfg.Color) && isequal(size(cfg.Color),[length(varargin),3])
    cfg.ColVal = cfg.Color;
else
    error('Unable to assign cfg.Color')
end

hold on;

ylims = get(gca,'Ylim');
lgnd = cell(size(varargin));
h = zeros(1,length(varargin)+1);

for iArg = 1:length(varargin)
    temp = varargin{iArg};
    if cfg.showLegend
        name = inputname(iArg+1); if isempty(name); name = ['Input ',num2str(iArg)]; end; lgnd{iArg} = name;
    end
    
    % check if Y is an iv data type
    if isfield(temp,'tstart')
       temp2 = [temp.tstart; temp.tend];
       x = [temp2 temp2 nan(size(temp2))]; 
       x = reshape(x',length(temp2)*3,1);
       
       y = repmat([ylims nan],length(temp2),1);
       y = reshape(y',length(temp2)*3,1);
       h(iArg) = plot(x,y,'Color',cfg.ColVal(iArg,:),'LineWidth',cfg.LineWidth);
       if cfg.patch
           patchY = []; patchX = []; % clear these from previous loop
           patchY = [ylims(1); ylims(2); ylims(2); ylims(1)];
           patchY = repmat(patchY,1,length(temp.tstart));
           patchX(1,:) = temp.tstart; patchX(2,:) = temp.tstart;
           patchX(3,:) = temp.tend; patchX(4,:) = temp.tend;
           
           patch(patchX,patchY,cfg.ColVal(iArg,:),'EdgeColor','none','FaceAlpha',0.1)
          
       end
    else
       h(iArg) = plot(temp.tvec,temp.data,'Color',cfg.ColVal(iArg,:),'LineWidth',cfg.LineWidth);
        
    end
end

% now do legend
if cfg.showLegend; h(iArg+1) = legend(h(1:iArg),lgnd); set(h(iArg+1),'Location',cfg.Location); end

end
