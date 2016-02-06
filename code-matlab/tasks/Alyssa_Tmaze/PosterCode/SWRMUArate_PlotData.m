%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%       Plot rate of SWR-associated MUA during session epochs         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Tmaze data
%
% run the script SWRMUArate_GenData before running this script.
%
%
%
%
% aacarey Sept 2015

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% do you want to save the figure?
cfg.writeFiles = 1; % 1 yes, 0 no

% where to save it and what to call it
cfg.output_fd = 'E:\Documents\TmazePaper\visuals';
cfg.output_fn = 'SWRMUArate';

% where to find the data to plot and what it's called
cfg.input_fd = 'E:\Documents\TmazePaper\data';
cfg.input_fn = 'SWRMUArate'; % do not specify extension, script knows

cfg.colormode = 'single'; 
% colormodes:
%    'inventory1': indiv rats are inventory colors; combined is grey
%    'inventory2': indiv rats are inventory; combined is red and blue
%    'rb': everything in red and blue
%    'grey': everything in grey

cfg.FontSize = 20;

%% Load the data

iWasHere = pwd; % remember where we started

cd(cfg.input_fd)
load(FindFile([cfg.input_fn,'.mat'])); % loads struct called 'precRate' for "precandidate rate"

% get colors to plot with
colors = TmazeColors(cfg.colormode); % colors is truct with food day and water day 
% colors for each rat. In some cases they might all be the same, but others 
% they are all different. this allows flexible plotting of different colors, if wanted

%% do the thing

rats = {'R042','R044','R050','R064','all'};

figure; hold on

for iRat = 1:length(rats)
   
   col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).f}; 
   d = [precRate.(rats{iRat}).pre precRate.(rats{iRat}).taskrest precRate.(rats{iRat}).post];

   location = [iRat iRat+length(rats)+2 iRat+2*length(rats)+4];
   
   for iBar = 1:length(d)
       h(iRat) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','k','LineWidth',2);
       hold on
   end
   
end

legtxt{1} = 'R042'; legtxt{2} = 'R044'; legtxt{3} = 'R050'; legtxt{4} = 'R064'; legtxt{5} = 'all rats'; 
legend(h(1:5),legtxt,'Location','north')

set(gca,'XTick',[],'FontSize',cfg.FontSize)
xlabel('     PRE               TASK REST             POST    ')
ylabel('Average rate of SWR events (n/min)')


% save it
cd(cfg.output_fd)
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-dpdf','-r300',cfg.output_fn);
print(gcf,'-depsc','-r300',cfg.output_fn);
cd(iWasHere)