% ALL_Plot_DeqSeq
%
% plot output from ALL_Collect_DecSeq

% some appearance things
cfg.colormode = 'inventory3';
FontSize = 8;
cfg.input_fd = 'C:\temp'; %'D:\projects\AlyssaTmaze\resultsFiles';
cfg.output_fd = 'C:\temp\viz'; %'D:\projects\AlyssaTmaze\resultsFiles\viz';
cfg.showAllRatsText = 1; % do you want to show the text "all rats" on the combined data figures?
cfg.writeOutput = 0;
cfg.input_prefix = 'R0_'; % which files to load? assumes filenames are *DecSeq_$task-phase_all_out.mat
cfg.outbasefn = 'R0_'; % base filename for figure output
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

%% load the data
cd(cfg.input_fd)

all = load(cat(2,cfg.input_prefix,'DecSeq_all_all_out'));
pre = load(cat(2,cfg.input_prefix,'DecSeq_prerecord_all_out'));
task = load(cat(2,cfg.input_prefix,'DecSeq_taskrest_all_out'));
post = load(cat(2,cfg.input_prefix,'DecSeq_postrecord_all_out'));

%% F0a: session-by-session averaged for all rats and task phases
%cfg.output_fn = cat(2,cfg.outbasefn,'sessionProps_overall');

% now collect data
this_data = all.data.all.ALL_sig_seq;

sess_label = repmat(1:6,[1 4]);
rat_label = cat(2,ones(1,6),2*ones(1,6),3*ones(1,6),4*ones(1,6));

sess_idx = sess_label(this_data.sess); % gives within-animal day idxs (1-6)
rat_idx = rat_label(this_data.sess);

% normalized (N) input data is redundant for left and right, so remove
sess_idx = sess_idx(1:2:end); rat_idx = rat_idx(1:2:end);
seq = this_data.countN(1:2:end);
behav = this_data.allTrialsN(1:2:end);
choice = this_data.choiceN(1:2:end);
type = this_data.type(1:2:end);

% split sessions into 'biased' and 'unbiased' behavior
behav_bias = max(cat(2,behav,1-behav),[],2); % invert normalized behavior below 0.5 to above 0.5
[behav_biasS,behav_bias_idx] = sort(behav_bias);
%behav_biasS(9:10); this suggests a threshold of 0.7 to divide
%biased/unbiased sessions
bias_thr = 0.7;
biased_idx = behav_bias_idx(behav_biasS > 0.7);
unbiased_idx = behav_bias_idx(behav_biasS <= 0.7);

% now need to find what the next session idxs are for these (may not
% exist)
what = {'biased_idx','unbiased_idx'};
for iW = 1:length(what)
    eval(cat(2,'these_idxs = ',what{iW}));
    temp_idx_out = {}; % this will hold idx into original data
    for iI = 1:length(these_idxs)
        this_idx = these_idxs(iI);
        this_sess = sess_idx(this_idx);
        this_rat = rat_idx(this_idx);
        
        %if this_sess == 6 | (this_rat == 2 & this_sess == 4) % last session for this rat, no next available
        %    temp_idx_out(iI) = [];
        %    continue;
        %end
        temp_idx_out{iI} = find(rat_idx == this_rat & sess_idx == this_sess + 1);
        
    end
    
    eval(cat(2,what{iW},'_out = cell2mat(temp_idx_out)'));
end

% assemble for plotting
figure;
what = {'biased_idx_out','unbiased_idx_out'};
for iW = 1:length(what)
    
    eval(cat(2,'this_idx = ',what{iW}));
    this_idx = intersect(this_idx,find(type == 1));
    food_left = nanmean(seq(this_idx)); food_right = nanmean(1-seq(this_idx));
    
    eval(cat(2,'this_idx = ',what{iW}));
    this_idx = intersect(this_idx,find(type == 2));
    water_left = nanmean(seq(this_idx)); water_right = nanmean(1-seq(this_idx));
    
    bar_idx = [1 2 4 5];
    subplot(2,3,iW);
    bar(bar_idx,[food_left food_right water_left water_right]);
    title(what{iW});
end



%% plot something



%% get averages for each session, and plot each data point individually
for iSess = 1:6
   sess_out(iSess) = iSess;
   seq_out(iSess) = nanmean(seq(sess_idx == iSess));
   behav_out(iSess) = nanmean(behav(sess_idx == iSess));
   
   plot(iSess-offs,seq(sess_idx == iSess),'.','MarkerSize',5,'Color',[0 0.7 0]);
   plot(iSess+offs,behav(sess_idx == iSess),'.','MarkerSize',5,'Color',[0 0 0]);
end

plot([0.5 6.5],[0.5 0.5],'LineStyle','--','LineWidth',1,'Color',[1 1 1]);
hold on;
set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
box off;

yyaxis right; set(gca,'YColor',[0 0.7 0],'LineWidth',1,'FontSize',fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'YLim',[0 1]);
ylabel('SWR sequence p(food)');


plot(sess_out+offs,behav_out,'k','LineStyle','-');
plot(sess_out+offs,behav_out,'.k','MarkerSize',20);

plot(sess_out-offs,seq_out,'Color',[0 0.7 0],'LineStyle','-');
plot(sess_out-offs,seq_out,'.','MarkerSize',20,'Color',[0 0.7 0]);

% make nice
set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
box off;

yyaxis left;
ylabel('choice p(food)');

if cfg.writeOutput
    maximize; drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end

% print some stats to go with this figure
tbl = table(categorical(rat_idx)',categorical(sess_idx)',behav,categorical(type),seq,choice,'VariableNames',{'Subject','Session','Behav','MotType','SeqContent','Choice'});
lme = fitglme(tbl,'SeqContent ~ MotType + (1|Subject)','Link','logit'); % subject-specific intercepts don't actually help!
lme2 = fitglme(tbl,'SeqContent ~ 1 + (1|Subject)','Link','logit');
compare(lme2,lme)

[r,p] = corrcoef(tbl.Choice,tbl.SeqContent)

%% F0b: session-by-session data for individual rats and task phases
cfg.output_fn = cat(2,cfg.outbasefn,'sessionProps_byPhase');

figure;

rats = {'R042','R044','R050','R064'};
what = {'pre','task','post'};

for iRat = 1:length(rats)
    for iW = 1:length(what)
        
        this_data = eval(cat(2,what{iW},'.data.',rats{iRat},'.ALL_sig_seq'));

        subplot(4,3,(iRat-1)*3+iW);
        set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
        box off;
        
        yyaxis right; set(gca,'YColor',[0 0.7 0],'LineWidth',1,'FontSize',fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'YLim',[0 1]);
      
        for iC = 1:6
            
            h(iC) = rectangle('Position',[iC-0.5 0 1 1]);
            if mod(iC,2) % odd
                if iRat >= 3
                    set(h(iC),'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                else
                    set(h(iC),'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                end
            else % even
                if iRat >= 3
                    set(h(iC),'FaceColor',[0.7 0.7 1],'EdgeColor',[0.7 0.7 1]);
                else
                    set(h(iC),'FaceColor',[1 0.7 0.7],'EdgeColor',[1 0.7 0.7]);
                end
            end
            
        end
        hold on;
        plot([0.5 6.5],[0.5 0.5],'LineStyle','--','LineWidth',1,'Color',[1 1 1]);
        set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
        
        plot(this_data.sess(1:2:end),this_data.allTrialsN(1:2:end),'k','LineStyle','-');
        plot(this_data.sess(1:2:end),this_data.allTrialsN(1:2:end),'.k','MarkerSize',20);
        
        plot(this_data.sess(1:2:end),this_data.countN(1:2:end),'Color',[0 0.7 0],'LineStyle','-');
        hold on;
        plot(this_data.sess(1:2:end),this_data.countN(1:2:end),'.','MarkerSize',20,'Color',[0 0.7 0]);
               
        % make nice
        set(gca,'XTick',1:6,'LineWidth',1,'TickDir','out','FontSize',fs,'YTick',0:0.25:1,'YTickLabel',{'0','','0.5','','1'},'XLim',[0.5 6.5],'YLim',[0 1]);
        box off;
          
        if iRat == 1
            title(what{iW});
        end
        
    end
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end
%% F1: raw sequence counts
cfg.output_fn = cat(2,cfg.outbasefn,'counts');

cfg.ylim = [450 250]; cfg.ylimtick = [150 125]; % ALL - overall and single lims
%cfg.ylim = [300 150]; cfg.ylimtick = [75 37.5]; % ALL - overall and single lims

ylab = {'Number of significant'; 'sequences'};
ylimsall = [0 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [0 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE
d = [pre.data.(rats{iRat}).food_left pre.data.(rats{iRat}).food_right pre.data.(rats{iRat}).water_left pre.data.(rats{iRat}).water_right];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',yticksall)
xlabel('  food                  water','FontSize',FontSize)
ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
title('PRERECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[3 4 9 10]) % for all rats TASKREST
d = [task.data.(rats{iRat}).food_left task.data.(rats{iRat}).food_right task.data.(rats{iRat}).water_left task.data.(rats{iRat}).water_right];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('TASK')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[5 6 11 12]) % for all rats POST
d = [post.data.(rats{iRat}).food_left post.data.(rats{iRat}).food_right post.data.(rats{iRat}).water_left post.data.(rats{iRat}).water_right];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('POSTRECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.65; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    d = [pre.data.(rats{iRat}).food_left pre.data.(rats{iRat}).food_right pre.data.(rats{iRat}).water_left pre.data.(rats{iRat}).water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    if position(1) == 13 || position(1) == 19
        yticks = ytickssing;
        ylabel(ylab,'FontSize',FontSize)
    else
        yticks = [];
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',yticks)
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(2)) % for indiv rat TASK
    d = [task.data.(rats{iRat}).food_left task.data.(rats{iRat}).food_right task.data.(rats{iRat}).water_left task.data.(rats{iRat}).water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(3)) % for indiv rat POST
    d = [post.data.(rats{iRat}).food_left post.data.(rats{iRat}).food_right post.data.(rats{iRat}).water_left post.data.(rats{iRat}).water_right];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end

%% F2: proportional sequence counts
cfg.ylim = [1 1]; cfg.ylimtick = [0.25 0.25]; % overall and single lims

colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_props');

ylab = {'Proportion of significant'; 'sequences'};
ylimsall = [0 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [0 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE
d = [pre.data.(rats{iRat}).food_leftN pre.data.(rats{iRat}).food_rightN pre.data.(rats{iRat}).water_leftN pre.data.(rats{iRat}).water_rightN];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',yticksall)
xlabel('  food                  water','FontSize',FontSize)
ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
title('PRERECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[3 4 9 10]) % for all rats TASKREST
d = [task.data.(rats{iRat}).food_leftN task.data.(rats{iRat}).food_rightN task.data.(rats{iRat}).water_leftN task.data.(rats{iRat}).water_rightN];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('TASK')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[5 6 11 12]) % for all rats POST
d = [post.data.(rats{iRat}).food_leftN post.data.(rats{iRat}).food_rightN post.data.(rats{iRat}).water_leftN post.data.(rats{iRat}).water_rightN];
for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end

set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',[])
xlabel('  food                  water','FontSize',FontSize)
title('POSTRECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.65; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    d = [pre.data.(rats{iRat}).food_leftN pre.data.(rats{iRat}).food_rightN pre.data.(rats{iRat}).water_leftN pre.data.(rats{iRat}).water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    if position(1) == 13 || position(1) == 19
        yticks = ytickssing;
        ylabel(ylab,'FontSize',FontSize)
    else
        yticks = [];
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',yticks)
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(2)) % for indiv rat TASK
    d = [task.data.(rats{iRat}).food_leftN task.data.(rats{iRat}).food_rightN task.data.(rats{iRat}).water_leftN task.data.(rats{iRat}).water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(3)) % for indiv rat POST
    d = [post.data.(rats{iRat}).food_leftN post.data.(rats{iRat}).food_rightN post.data.(rats{iRat}).water_leftN post.data.(rats{iRat}).water_rightN];
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimssing,'YTick',[])
    xlabel('  food   water','FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end

%% F3: relative differences
colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_props_rel');

cfg.xlim = [0.5 0.8]; cfg.xlimtick = [0.5 0.8]; % overall and single lims

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).w}; % {foodColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE

food_diff = pre.data.(rats{iRat}).food_rightN-pre.data.(rats{iRat}).food_leftN;
water_diff = pre.data.(rats{iRat}).water_rightN-pre.data.(rats{iRat}).water_leftN;

h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
hold on
h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(1) cfg.xlim(1)],'XTick',-cfg.xlim(1):cfg.xlimtick(1)/2:cfg.xlim(1), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(1)),-cfg.xlim(1)/2,0,cfg.xlim(1)/2,sprintf('W %.1f',cfg.xlim(1))}, ...
    'YTick',1:2,'YTickLabel',{'WaterR','FoodR'},'FontSize',FontSize,'YLim',[0.5 2.5]);

title('PRERECORD')
box off
set(gca,'Layer','top')

if cfg.showAllRatsText
txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[3 4 9 10]) % for all rats TASKREST
food_diff = task.data.(rats{iRat}).food_rightN-task.data.(rats{iRat}).food_leftN;
water_diff = task.data.(rats{iRat}).water_rightN-task.data.(rats{iRat}).water_leftN;

h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
hold on
h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(1) cfg.xlim(1)],'XTick',-cfg.xlim(1):cfg.xlimtick(1)/2:cfg.xlim(1), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(1)),-cfg.xlim(1)/2,0,cfg.xlim(1)/2,sprintf('W %.1f',cfg.xlim(1))}, ...
    'YTick',1:2,'YTickLabel',{'WaterR','RoodR'},'FontSize',FontSize,'YLim',[0.5 2.5]);

title('TASK')
box off
set(gca,'Layer','top')

if cfg.showAllRatsText
txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

subplot(4,6,[5 6 11 12]) % for all rats POST

food_diff = post.data.(rats{iRat}).food_rightN-post.data.(rats{iRat}).food_leftN;
water_diff = post.data.(rats{iRat}).water_rightN-post.data.(rats{iRat}).water_leftN;

h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
hold on
h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(1) cfg.xlim(1)],'XTick',-cfg.xlim(1):cfg.xlimtick(1)/2:cfg.xlim(1), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(1)),-cfg.xlim(1)/2,0,cfg.xlim(1)/2,sprintf('W %.1f',cfg.xlim(1))}, ...
    'YTick',1:2,'YTickLabel',{'WaterR','FoodR'},'FontSize',FontSize,'YLim',[0.5 2.5]);

title('POST')
box off
set(gca,'Layer','top')

if cfg.showAllRatsText
txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.2; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).w}; % {foodColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    food_diff = pre.data.(rats{iRat}).food_rightN-pre.data.(rats{iRat}).food_leftN;
    water_diff = pre.data.(rats{iRat}).water_rightN-pre.data.(rats{iRat}).water_leftN;
    
    h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
    hold on
    h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(2) cfg.xlim(2)],'XTick',-cfg.xlim(2):cfg.xlimtick(2)/2:cfg.xlim(2), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(2)),'',0,'',sprintf('W %.1f',cfg.xlim(2))}, ...
    'YTick',1:2,'YTickLabel',{'Wr','Fr'},'FontSize',FontSize,'YLim',[0.5 2.5]);
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(2)) % for indiv rat TASK
    
    food_diff = task.data.(rats{iRat}).food_rightN-task.data.(rats{iRat}).food_leftN;
    water_diff = task.data.(rats{iRat}).water_rightN-task.data.(rats{iRat}).water_leftN;
    
    h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
    hold on
    h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(2) cfg.xlim(2)],'XTick',-cfg.xlim(2):cfg.xlimtick(2)/2:cfg.xlim(2), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(2)),'',0,'',sprintf('W %.1f',cfg.xlim(2))}, ...
    'YTick',1:2,'YTickLabel',{'Wr','Fr'},'FontSize',FontSize,'YLim',[0.5 2.5]);
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
    subplot(4,6,position(3)) % for indiv rat POST
    
    food_diff = post.data.(rats{iRat}).food_rightN-post.data.(rats{iRat}).food_leftN;
    water_diff = post.data.(rats{iRat}).water_rightN-post.data.(rats{iRat}).water_leftN;
    
    h(1) = barh(1,nanmean(water_diff)); set(h(1),'FaceColor',col{2},'EdgeColor','none');
    hold on
    h(2) = barh(2,nanmean(food_diff)); set(h(2),'FaceColor',col{1},'EdgeColor','none');
set(gca,'XLim',[-cfg.xlim(2) cfg.xlim(2)],'XTick',-cfg.xlim(2):cfg.xlimtick(2)/2:cfg.xlim(2), ...
    'XTickLabel',{sprintf('%.1f F',-cfg.xlim(2)),'',0,'',sprintf('W %.1f',cfg.xlim(2))}, ...
    'YTick',1:2,'YTickLabel',{'Wr','Fr'},'FontSize',FontSize,'YLim',[0.5 2.5]);
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end


%% F4: decoding error
cfg.ylim = [20 30]; cfg.ylimtick = [5 10]; % overall and single lims

colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_decAcc');

pre = load(cat(2,cfg.input_prefix,'DecSeq_prerecord_all_out')); % why reload?

ylab = {'decoding error'};

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE

this_data = pre.data.(rats{iRat}).ALL_sig_seq;

food_left_idx = find(this_data.arm == 1 & this_data.type == 1);
d(1) = nanmean(this_data.decErr(food_left_idx));
food_right_idx = find(this_data.arm == 2 & this_data.type == 1);
d(2) = nanmean(this_data.decErr(food_right_idx));
water_left_idx = find(this_data.arm == 1 & this_data.type == 2);
d(3) = nanmean(this_data.decErr(water_left_idx));
water_right_idx = find(this_data.arm == 2 & this_data.type == 2);
d(4) = nanmean(this_data.decErr(water_right_idx));

for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',[0 cfg.ylim(1)],'YTick',0:cfg.ylimtick(1):cfg.ylim(1));
xlabel('  food                  water','FontSize',FontSize)
ylabel(ylab,'FontSize',FontSize)
title('DecAcc')
box off

set(gca,'Layer','top')
if cfg.showAllRatsText
    txt = 'all rats';
    text(0.15,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.2; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
      
    subplot(4,6,position(1)) % for indiv rat PRE
    
    
    this_data = pre.data.(rats{iRat}).ALL_sig_seq;
    
    food_left_idx = find(this_data.arm == 1 & this_data.type == 1);
    d(1) = nanmean(this_data.decErr(food_left_idx));
    food_right_idx = find(this_data.arm == 2 & this_data.type == 1);
    d(2) = nanmean(this_data.decErr(food_right_idx));
    water_left_idx = find(this_data.arm == 1 & this_data.type == 2);
    d(3) = nanmean(this_data.decErr(water_left_idx));
    water_right_idx = find(this_data.arm == 2 & this_data.type == 2);
    d(4) = nanmean(this_data.decErr(water_right_idx));
    
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',[0 cfg.ylim(2)],'YTick',0:cfg.ylimtick(2):cfg.ylim(2));
    xlabel('  food                  water','FontSize',FontSize)
    ylabel(ylab,'FontSize',FontSize)
    box off
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end

%% F5: arranged by first choice
cfg.ylim = [1 1]; cfg.ylimtick = [0.25 0.25]; % overall and single lims

colors = TmazeColors(cfg.colormode);

originalFolder = pwd;

cd(cfg.input_fd)
cfg.output_fn = cat(2,cfg.outbasefn,'_firstchoice');

ylab = {'Proportion of significant'; 'sequences'};
ylimsall = [0 cfg.ylim(1)];
yticksall = ylimsall(1):cfg.ylimtick(1):ylimsall(2);
ylimssing = [0 cfg.ylim(2)];
ytickssing = ylimssing(1):cfg.ylimtick(2):ylimssing(2);

location = [1 2 4 5]; % where to place the bar
xlims = [0 location(4)+1];

% get the accumulated data subplotted first 
 
rats = {'all'};
iRat = 1;

col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
figure; hold on

subplot(4,6,[1 2 7 8]) % for all rats PRE

this_data = pre.data.(rats{iRat}).ALL_sig_seq;

choiceleft_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 1);
choiceleft_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 1);

choiceright_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 2);
choiceright_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 2);


xbin = [1 2 4 5];
d = [nanmean(choiceleft_left) nanmean(choiceleft_right) nanmean(choiceright_left) nanmean(choiceright_right)];
cl = 'rrbb';

for iBar = 1:length(d)
    h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
    hold on
end
set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
set(gca,'YLim',ylimsall,'YTick',yticksall)
xlabel('  left                  right','FontSize',FontSize)
ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
title('PRERECORD')
box off
set(gca,'Layer','top')
if cfg.showAllRatsText
txt = 'all rats';
    text(0.58,0.9,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
end

% individual rat data
rats = {'R042','R044','R050','R064'}; % DO NOT CHANGE
%rats = {'R042','R050','R064'};

start = [13 14 19 20]; % R042,R044,R050,R064 start positions for PRE
xName = 0.2; yName = 0.85; % controls where to place rat name text

for iRat = 1:length(rats)
    
    position = [start(iRat) start(iRat)+2 start(iRat)+4];
    col = {colors.(rats{iRat}).f colors.(rats{iRat}).f colors.(rats{iRat}).w colors.(rats{iRat}).w}; % {foodColor foodColor waterColor waterColor}
    
    subplot(4,6,position(1)) % for indiv rat PRE
    
    this_data = pre.data.(rats{iRat}).ALL_sig_seq;
    
    choiceleft_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 1);
    choiceleft_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 1);
    
    choiceright_left = this_data.countN(this_data.arm == 1 & this_data.firstChoice == 2);
    choiceright_right = this_data.countN(this_data.arm == 2 & this_data.firstChoice == 2);
    
    
    xbin = [1 2 4 5];
    d = [nanmean(choiceleft_left) nanmean(choiceleft_right) nanmean(choiceright_left) nanmean(choiceright_right)];
    cl = 'rrbb';
    
    for iBar = 1:length(d)
        h(iBar) = bar(location(iBar),d(iBar),'FaceColor',col{iBar},'EdgeColor','none');
        hold on
    end
    set(gca,'XTick',location,'XTickLabel',{'L' 'R' 'L' 'R'},'FontSize',FontSize,'LineWidth',1,'XLim',xlims);
    set(gca,'YLim',ylimsall,'YTick',yticksall)
    xlabel('  left                  right','FontSize',FontSize)
    ylabel([ylab{1},' ',ylab{2}],'FontSize',FontSize)
    box off
    set(gca,'Layer','top')
    
    box off
    set(gca,'Layer','top')
    txt = rats{iRat};
    text(xName,yName,txt,'Units','normalized',...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','right',...
        'FontSize',FontSize)
    
end

% save thing
if cfg.writeOutput
    maximize; set(gcf,'PaperPositionMode','auto'); drawnow;
    cd(cfg.output_fd);
    print(gcf,'-painters','-dpng','-r300',[cfg.output_fn,'.png']);
    print(gcf,'-painters','-dpdf','-r300',[cfg.output_fn,'.pdf']);
    print(gcf,'-painters','-depsc','-r300',[cfg.output_fn,'.eps']);
    cd(originalFolder)
end