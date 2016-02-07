function ax = PolyRaster(cfg,S1,varargin)
%POLYRASTER Plot multiple spike rasters and link x axes for synchronous
%navigation.
%   POLYRASTER(cfg,S1,varargin)
%   POLYRASTER(cfg,S1,S2,S3, etc)
%
%   Calls on MultiRaster() in subplotted axes. See navigate() for keyboard
%   commands.
%
%   INPUTS
%      cfg: - single struct with fields controlling function behaviour, OR
%           - {1 x nS} cell array containing configs in the same order as your
%             inputs (these will be paired with the corresponding S during
%             plotting)
%       S1: TS struct containing spiketrains
%       S2 and so on: other spiketrains from the same session, perhaps with
%           a different order for comparison
%
%   OUTPUTS
%       ax: axes handles, ax(1) corresponds to S1, ax(2) to S2 and so on.
%           Use the handles to set figure titles and stuff
%           ex: set(get(ax(2),'title'),'string','Left spiketrains')
%
%   CONFIG OPTIONS
%       See MultiRaster
%
%    See also: MultiRaster, PlotSpikeRaster2, navigate 
%
% aacarey Jan 2016

nS = 1 + length(varargin);

% Make sure the user didn't mess up config numbers
if length(cfg)>1
    assert(length(cfg) == nS,'If cfg is an array, it must contain the same number of configs as there are S inputs')
    nConfig = 'multi';
else
    nConfig = 'single';
end

% Collect S data into an array to facilitate loop plotting
allS(1) = S1;
for iVarg = length(varargin):-1:1
    allS(iVarg+1) = varargin{iVarg};
end

% Loop plot
figure('Name',mfilename,'KeyPressFcn',@navigate)
for iS = nS:-1:1
    ax(iS) = subplot(nS,1,iS);
    switch nConfig
        case 'single'
            cfg.openNewFig = 0;
            MultiRaster(cfg,allS(iS));
        case 'multi'
            cfg{1,iS}.openNewFig = 0;
            MultiRaster(cfg{1,iS},allS(iS));
    end
end

% link the x axes for synchronous navigation
linkaxes(ax,'x')

end

