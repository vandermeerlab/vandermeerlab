%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%          Plot Prerecord, Task, Postrecord CoOccurrence Data         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For Tmaze project
%
% Plots Co-Occurrence data in a different way. Gets session average for
% CoOccurQ p numbers for food and water combined and then subtracts it from
% each p number. What is plotted is the resulting mean food numbers minus
% the resulting mean water numbers.
% A "processed" way of visualizing the co-occur data
% 

clearvars -except CFG

%% WHAT DO YOU WANT THIS SCRIPT TO DO?

% what to plot?
cfg.whichP = 'p4'; % 'p0' or 'p4' (single cell or zscore coactivity), can also do p1, p2, p3, p5 now

% some appearance things
cfg.colormode = 'inventory4';
FontSize = 8;

% where to save the figure?
cfg.input_fd = [pwd,'\data'];
cfg.output_fd = [pwd,'\visuals'];

cfg.whichITI = 'equalBehavior'; % 'equalBehavior' or 'all'

%% If master config exists, process cfg

if exist('CFG','var')
    cfg = ProcessConfig2(cfg,CFG.plot);
end

%% Set some things up

colors = TmazeColors(cfg.colormode); % get colours for plotting

originalFolder = pwd;

cd(cfg.input_fd)

% get prerecord co-occur data
data.pre = loadpop('coOccurrence_prerecord.mat','coocData');

switch cfg.whichITI
    case 'all'
        % get all ITI data        
        data.allITI = loadpop('coOccurrence_allITI.mat','coocData');
    case 'equalBehavior'
        % get equal behaviour ITI data
        data.equalBehaviorITI = loadpop('CoOccurrence_equalBehaviorITI.mat','coocData');
end

% get task data
data.task = loadpop('coOccurrence_task.mat','coocData');

% get postrecord co-occur data
data.post = loadpop('coOccurrence_postrecord.mat','coocData');


%% Collect data 

% dynamic field names I'll be using:
epochs = fieldnames(data); % 'pre','task','equalBehaviourITIT', 'post'
rats = TmazeRats;
restrictionType = {'food','water'};
arms = {'L','R'};

% Get session means
for iEpoch = 1:length(epochs)
    for iRat = 1:length(rats)
        disp(' '); disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(['~~  COLLECTING CO-OCCURRENCE DATA FOR ',rats{iRat},' ~~'])
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        
        sessionType.food = []; sessionType.water = [];
        
        %~~~~~~~~~~~~~~~~~~~~ iterate through sessions ~~~~~~~~~~~~~~~~~~~~~~~~
        for iSession = 1:length(data.(epochs{iEpoch}).(rats{iRat}))
            
            % get restriction type session numbers
            switch data.(epochs{iEpoch}).(rats{iRat})(iSession).restrictionType
                case 'food'
                    sessionType.food = [sessionType.food,iSession];
                case 'water'
                    sessionType.water = [sessionType.water,iSession];
            end
            
            % get mean proportions for left and right together
            sessionMean = nanmean([data.(epochs{iEpoch}).(rats{iRat})(iSession).L.ALLp.(cfg.whichP); data.(epochs{iEpoch}).(rats{iRat})(iSession).R.ALLp.(cfg.whichP)]);
            
            %~~~~~~~~~~~~~~~~~ iterate through arms ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            for iArm = 1:length(arms)
                p = data.(epochs{iEpoch}).(rats{iRat})(iSession).(arms{iArm}).ALLp.(cfg.whichP);
                ccDiffbySession.(epochs{iEpoch}).(rats{iRat})(iSession).(arms{iArm}).p = nanmean(p - sessionMean);
            end
            
            % SUBTRACT L-R
            ccDiffbySession.(epochs{iEpoch}).(rats{iRat})(iSession).diff = ccDiffbySession.(epochs{iEpoch}).(rats{iRat})(iSession).R.p - ccDiffbySession.(epochs{iEpoch}).(rats{iRat})(iSession).L.p;
        end
        
        % get contingency means (not sure if good word) 'foodL' 'foodR' 'waterL' 'waterR'
        
        for iRestr = 1:length(restrictionType)
            
            ccDiff.(epochs{iEpoch}).(rats{iRat}).(restrictionType{iRestr}) = nanmean([ccDiffbySession.(epochs{iEpoch}).(rats{iRat})(sessionType.(restrictionType{iRestr})(1)).diff,...
                ccDiffbySession.(epochs{iEpoch}).(rats{iRat})(sessionType.(restrictionType{iRestr})(2)).diff,...
                ccDiffbySession.(epochs{iEpoch}).(rats{iRat})(sessionType.(restrictionType{iRestr})(3)).diff]);
        end
    end % of rats
    
    % coolect for avg of rat data
    for iRestr = 1:length(restrictionType) % got lazy, hard-coded the rat fields... :(
            ccDiff.(epochs{iEpoch}).all.(restrictionType{iRestr}) =[ccDiff.(epochs{iEpoch}).R042.(restrictionType{iRestr})...
                ccDiff.(epochs{iEpoch}).R044.(restrictionType{iRestr})...
                ccDiff.(epochs{iEpoch}).R050.(restrictionType{iRestr})...
                ccDiff.(epochs{iEpoch}).R064.(restrictionType{iRestr})];
            
            % get SEM values
            ysemSubj.(epochs{iEpoch}).all.(restrictionType{iRestr}) = nanstd(ccDiff.(epochs{iEpoch}).all.(restrictionType{iRestr}))/sqrt(length(rats));
            ysemSess.(epochs{iEpoch}).all.(restrictionType{iRestr}) = nanstd(ccDiff.(epochs{iEpoch}).all.(restrictionType{iRestr}))/sqrt(length(rats)*3);
            ccDiff.(epochs{iEpoch}).all.(restrictionType{iRestr}) = nanmean(ccDiff.(epochs{iEpoch}).all.(restrictionType{iRestr}));
    end
end

%% PLOTTING SECTION

% general
tickBuffer = 0.01;

switch cfg.whichP
    case 'p0'
        xlab = {'Proportion of';'SWRs active'};
        cfg.output_fn = 'p0DiffMulti';
        
    case 'p1'
        xlab = {'Cell pair'; 'whatever this is'};
        cfg.output_fn = 'p1DiffMulti';
        
    case 'p2'
        xlab = {'Cell pair'; 'conditional probability'};
        cfg.output_fn = 'p2DiffMulti';
        
    case 'p3'
        xlab = {'Cell pair'; 'joint probability'};
        cfg.output_fn = 'p3DiffMulti';
        
    case 'p4'
        xlab = {'SWR coactivation'; 'Z-score'};
        cfg.output_fn = 'p4DiffMulti';
        
    case 'p5'
        xlab = {'Randomly shuffled'; 'stuff'};
        cfg.output_fn = 'p5DiffMulti';
        
end


location = [1 2]; % where to place the bar
xlims = [0 location(2)+1];
epochs = fieldnames(ccDiff);

% add the 'all' category
rats = ['all',rats];

% specific to "all rats" group

position_all = [1 2 9 10; 3 4 11 12; 5 6 13 14; 7 8 15 16];
start = [17 18 25 26];
jhgjhg = [0 2 4 6];

figure; hold on

for iEpoch = 1:length(epochs)
    for iRat = 1:length(rats)
        switch rats{iRat}
            case 'all'
                position = position_all(iEpoch,:);
                txt = 'all rats';
            otherwise
                position = start(iRat-1) + jhgjhg(iEpoch); % tracks which spot the subplot will occupy
                txt = rats{iRat};
        end
        
        col = {colors.(rats{iRat}).f colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
        
        subplot(4,8,position)
        d = [ccDiff.(epochs{iEpoch}).(rats{iRat}).food ccDiff.(epochs{iEpoch}).(rats{iRat}).water];
        for iBar = 1:length(d)
            bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
            hold on
        end
        
        if iRat == 1
%             errorbar(location-0.1,d,[ysemSubj.(epochs{iEpoch}).all.food ysemSubj.(epochs{iEpoch}).all.water],'k','LineStyle','none');
%             errorbar(location+0.1,d,[ysemSess.(epochs{iEpoch}).all.food ysemSess.(epochs{iEpoch}).all.water],'Color',[0.5 0.5 0.5],'LineStyle','none');
            errorbar(location,d,[ysemSess.(epochs{iEpoch}).all.food ysemSess.(epochs{iEpoch}).all.water],'Color','k','LineStyle','none');

            % store some info pertaining to the plots for formatting later
            h.(rats{iRat})(iEpoch) = gca;
            mindata.(epochs{iEpoch})(iRat) = min([d(1) - ysemSubj.(epochs{iEpoch}).all.food,d(1) - ysemSess.(epochs{iEpoch}).all.food, d(2) - ysemSubj.(epochs{iEpoch}).all.water, d(2) - ysemSess.(epochs{iEpoch}).all.water]); % for y limits assignment, get lowest value
            maxdata.(epochs{iEpoch})(iRat) = max([d(1) + ysemSubj.(epochs{iEpoch}).all.food,d(1) + ysemSess.(epochs{iEpoch}).all.food, d(2) + ysemSubj.(epochs{iEpoch}).all.water, d(2) + ysemSess.(epochs{iEpoch}).all.water]); % for y limits, get highest value
            
            ylabel('    L                        R   ','FontSize',FontSize)
        else
            % store some info pertaining to the plots for formatting later
            h.(rats{iRat})(iEpoch) = gca;
            mindata.(epochs{iEpoch})(iRat) = min(d); % for y limits assignment, get lowest value
            maxdata.(epochs{iEpoch})(iRat) = max(d); % for y limits, get highest value
            
            ylabel(' L         R','FontSize',FontSize)
        end
        if isequal(position,[1 2 9 10]) || position(1) == 17 || position(1) == 25
            xlabel(xlab,'FontSize',FontSize)
        else
            yticks = [];
        end
        box off
        view(270,90); % rotate the axes
    
        if iRat == 1
            switch epochs{iEpoch}
                case 'pre'
                    title('PRERECORD','FontSize',FontSize)
                case 'task'
                    title('TASK','FontSize',FontSize)
                case 'equalBehaviorITI'
                    title('EQUAL BEHAVIOR ITI','FontSize',FontSize)
                case 'allITI'
                    title('ALL INTERTRIAL','FontSize',FontSize)
                case 'post'
                    title('POSTRECORD','FontSize',FontSize)
            end
        end
        
        if iRat == 1
            text(0.58,0.9,txt,'Units','normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','right',...
                'FontSize',FontSize,'BackgroundColor','w')
        else
            text(0.65,0.9,txt,'Units','normalized',...
                'VerticalAlignment','bottom',...
                'HorizontalAlignment','right',...
                'FontSize',FontSize,'BackgroundColor','w')
        end
    end % of rats
end % of epochs

for iRat = 1:length(rats)
    set(h.(rats{iRat}),'XTick',location,'XTickLabel',{'F' 'W'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims,'Layer','top');
    
    for iEpoch = 1:length(epochs)
        
        if iRat ~= 1
            
            extreme = max([abs(mindata.(epochs{iEpoch})) abs(maxdata.(epochs{iEpoch}))]);
            extreme = extreme + extreme*0.1;
            if ~isnan(extreme)
            set(h.(rats{iRat})(iEpoch), 'YLim', [-extreme extreme],'YDir','reverse')
            end
            
        else
            
            %extreme = max([abs(mindata.(epochs{iEpoch})(1)) %abs(maxdata.(epochs{iEpoch})(1))]); % each has unique ylims
            extreme = 1.05; % hardcode ylims, all same
            extreme = extreme + extreme*0.1;
            set(h.(rats{iRat})(iEpoch), 'YLim', [-extreme extreme],'YDir','reverse')
            
        end
        
        % make all yticks positive
%         currentTicks = get(h.(rats{iRat})(iEpoch),'ytick');
%         newTicks = num2str(abs(currentTicks.'));
%         set(h.(rats{iRat})(iEpoch),'YTick',currentTicks,'YTickLabel',newTicks)
        
    end % of epochs
end % of rats
hold off

%% If master config exists, save config history

if exist('CFG','var')
    CFG = History(CFG,mfilename,cfg);
end

%% save it
cd(cfg.output_fd)
print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
print(gcf,'-depsc','-r300',cfg.output_fn);
print(gcf,'-dpdf','-r300',cfg.output_fn);

cd(originalFolder)

disp('~~~~~~~ END OF PLOT RUN ~~~~~~~')