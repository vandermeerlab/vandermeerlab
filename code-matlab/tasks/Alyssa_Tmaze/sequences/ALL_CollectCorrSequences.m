% batch script to collect and plot co-occurrence data

%%
cfg = []; cfg.rats = {'R042'};
cfg.requireCandidates = 1;
fd = getTmazeDataPath(cfg);

cfg.indivSessionPlot = 1;
cfg.indivSequencePlots = 1;

cfg.prefix = 'CS_ALL_'; % which files to load
cfg.p = 0.99; % p-level for which events to count
cfg.sessions = {'food','water'};
cfg.arms = {'left','right'}; % needs to match order of out.score1[] and out.score3[] (1 is left, 2 right arm)
cfg.writeFiles = 1;
cfg.eventsToProcess = 'all'; % 'postrecord', 'all'

%% init vars to place data into

ALL_sig_seq.count = []; % the actual observed activation probabilities
ALL_sig_seq.arm = []; % left (1), right (2)
ALL_sig_seq.type = []; % restriction type (food/water)
ALL_sig_seq.sess = []; % session ID in case we want to restrict later

%% collect
for iFD = 1:length(fd)
    
    close all;
    
   cd(fd{iFD}); 
   LoadExpKeys;
   
   %cd('files');
   cd([fd{iFD},'\files'])
   [~,sessionID,~] = fileparts(fd{iFD}); % 'RXXX-201X-XX-XX'
   this_file = FindFiles([cfg.prefix,sessionID,'-CorrScores.mat']);
   
   if isempty(this_file)
      fprintf('Session %s: no CorrScores file found, skipping...\n',fd{iFD}); 
      continue;
   end
   
   load(this_file{1}); cd ..
   
   this_session_type = find(strcmp(ExpKeys.RestrictionType,cfg.sessions));
   
   % simply count the number of significant sequences
   for iArm = 1:length(cfg.arms)
       
       % optionally, restrict to specific task epoch
       
       keep_idx{iArm} = find(out.score1(iArm).WIN_rho_perc > cfg.p & out.score3(iArm).WIN_rho_perc > cfg.p);
   
       switch cfg.eventsToProcess
           case 'prerecord'
               t = [out.score1(iArm).WIN_iv(:).tstart];
               t = t(keep_idx{iArm});
               keep_idx{iArm} = keep_idx{iArm}(t < ExpKeys.prerecord(2));
           case 'postrecord'
               t = [out.score1(iArm).WIN_iv(:).tstart];
               t = t(keep_idx{iArm});
               keep_idx{iArm} = keep_idx{iArm}(t > ExpKeys.postrecord(1));
       end
                 
       ALL_sig_seq.count = cat(1,ALL_sig_seq.count,length(keep_idx{iArm}));
       ALL_sig_seq.arm = cat(1,ALL_sig_seq.arm,iArm);
       ALL_sig_seq.type = cat(1,ALL_sig_seq.type,this_session_type);
       ALL_sig_seq.sess = cat(1,ALL_sig_seq.sess,iFD);
       
   end
   
   % output something about this session
   fprintf('%s (%s): p0 L %d, R %d\n', ...
       ExpKeys.goodSWR{1}(1:15),cfg.sessions{this_session_type},length(keep_idx{1}),length(keep_idx{2}));
   
   % if specified, 
   if cfg.indivSessionPlot
       
       load(FindFile('*metadata.mat'));
       

       figure;
       
       % behavior
       a1 = subplot(311);
       plotcfg = []; plotcfg.ax = a1;
       PlotSessionBehavior(plotcfg,metadata,ExpKeys); 
       xl = get(gca,'XLim');
       
       % number of cells in each event
       spk = GetPC([]); % gets spk.S_pc(1) and (2) for left and right cells
       
       temp_tstart = [out.score1(1).WIN_iv.tstart];
       temp_tend = [out.score1(1).WIN_iv.tend];
       evtL = iv(temp_tstart,temp_tend);
       
       cfg_n = []; cfg_n.label = 'NActiveCells_left';
       evtL = AddNActiveCellsIV(cfg_n,evtL,spk.S_pc(1)); % note this has binning at 5ms by default
       
       temp_tstart = [out.score1(2).WIN_iv.tstart];
       temp_tend = [out.score1(2).WIN_iv.tend];
       evtR = iv(temp_tstart,temp_tend);
       
       cfg_n = []; cfg_n.label = 'NActiveCells_left';
       evtR = AddNActiveCellsIV(cfg_n,evtR,spk.S_pc(2)); % note this has binning at 5ms by default
       
       % sequences - by nCells
       a2 = subplot(312);
       
       all_tvec = [out.score1(1).WIN_iv.tstart]; % times of all sequences
       
       left_IDperc = out.score1(1).WIN_rho_perc; left_Tperc = out.score3(1).WIN_rho_perc;
       right_IDperc = out.score1(2).WIN_rho_perc; right_Tperc = out.score3(2).WIN_rho_perc;
       
       l95_idx = find(left_IDperc > 0.95 & left_Tperc > 0.95);
       l99_idx = find(left_IDperc > 0.99 & left_Tperc > 0.99);
       
       %plot(all_tvec(l95_idx),left_IDperc(l95_idx).*left_Tperc(l95_idx),'r.');
       %hold on;
       %plot(all_tvec(l99_idx),left_IDperc(l99_idx).*left_Tperc(l99_idx),'ro');
       
       plot(all_tvec(l95_idx),evtL.usr(1).data(l95_idx),'.','Color',[1 0.5 0.5],'MarkerSize',5);
       hold on;
       plot(all_tvec(l99_idx),evtL.usr(1).data(l99_idx),'.r','MarkerSize',20);
       plot(all_tvec(l99_idx),evtL.usr(1).data(l99_idx),'or','MarkerSize',10);
       
       r95_idx = find(right_IDperc > 0.95 & right_Tperc > 0.95);
       r99_idx = find(right_IDperc > 0.99 & right_Tperc > 0.99);
       
       %plot(all_tvec(r95_idx),1-(right_IDperc(r95_idx).*right_Tperc(r95_idx)),'b.');
       %plot(all_tvec(r99_idx),1-(right_IDperc(r99_idx).*right_Tperc(r99_idx)),'bo');
       
       plot(all_tvec(r95_idx),-evtR.usr(1).data(r95_idx),'.','Color',[0.5 0.5 1],'MarkerSize',5);
       plot(all_tvec(r99_idx),-evtR.usr(1).data(r99_idx),'.b','MarkerSize',20);
       plot(all_tvec(r99_idx),-evtR.usr(1).data(r99_idx),'bo','MarkerSize',10);
       
       set(gca,'XLim',xl,'YTick',-15:5:15,'YTickLabel',{'15','','',0,'','','15'},'YLim',[-15 15],'FontSize',14);
       ylabel('< R |  nCells  | L >'); box off;
       
       plot(xl,[5 5],'r:'); plot(xl,[-5 -5],'b:');
       
       if cfg.writeFiles
            cd([fd{iFD},'\files']) %cd('files');
           base_fn = ExpKeys.goodSWR{1}(1:15);
           
           fn = cat(2,base_fn,'_indivSessionCorr.png');
           print(gcf,'-dpng','-r300',fn);
           fn = cat(2,base_fn,'_indivSessionCorr.ai');
           print(gcf,'-dill',fn);
           
           cd ..
       end
       
   end
   
   %
   if cfg.indivSequencePlots
   
       % get pc_left and pc_right! should be a function...
       spk = GetPC([]);
       
       % L
       l99_iv = out.score1(1).WIN_iv(l99_idx);
       
       plot_cfg = []; plot_cfg.evt.tstart = [l99_iv.tstart]; plot_cfg.evt.tend = [l99_iv.tend];
       plot_cfg.LineWidth = 2; plot_cfg.axislabel = 'off';
       MultiRasterTwin(plot_cfg,spk.S_pc(1),spk.S_pc(2));
       
       for iEvt = 1:length(plot_cfg.evt.tstart)
           c = nanmean([plot_cfg.evt.tstart(iEvt) plot_cfg.evt.tend(iEvt)]);
           set(gca,'XLim',[c-0.5 c+0.5]);
           drawnow;
           
           if cfg.writeFiles
                cd([fd{iFD},'\files']) %cd('files');
               
               fn = cat(2,'IndivSeqL_',num2str(iEvt),'.png');
               print(gcf,'-dpng','-r300',fn);
               fn = cat(2,'IndivSeqL_',num2str(iEvt),'.ai');
               print(gcf,'-dill',fn);
               
               cd ..
           end
           
       end % over events (left)
       
       r99_iv = out.score1(2).WIN_iv(r99_idx);
       
       plot_cfg = []; plot_cfg.evt.tstart = [r99_iv.tstart]; plot_cfg.evt.tend = [r99_iv.tend];
       plot_cfg.LineWidth = 2; plot_cfg.axislabel = 'off';
       MultiRasterTwin(plot_cfg,spk.S_pc(1),spk.S_pc(2));
       
       for iEvt = 1:length(plot_cfg.evt.tstart)
           c = nanmean([plot_cfg.evt.tstart(iEvt) plot_cfg.evt.tend(iEvt)]);
           set(gca,'XLim',[c-0.5 c+0.5]);
           drawnow;
           
           if cfg.writeFiles
                cd([fd{iFD},'\files']) %cd('files');
               
               fn = cat(2,'IndivSeqR_',num2str(iEvt),'.png');
               print(gcf,'-dpng','-r300',fn);
               fn = cat(2,'IndivSeqR_',num2str(iEvt),'.ai');
               print(gcf,'-dill',fn);
               
               cd ..
           end
           
       end % over events (right)
       
       
       
   end
   
end

%% summary plot
figure; 

food_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 1);
food_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 1);

water_left = ALL_sig_seq.count(ALL_sig_seq.arm == 1 & ALL_sig_seq.type == 2);
water_right = ALL_sig_seq.count(ALL_sig_seq.arm == 2 & ALL_sig_seq.type == 2);

subplot(221);

xbin = [1 2 4 5];
d = [nansum(food_left) nansum(food_right) nansum(water_left) nansum(water_right)];
cl = 'rrbb';

for iB = 1:length(xbin)
   
    h(iB) = bar(xbin(iB),d(iB)); set(h(iB),'FaceColor',cl(iB),'EdgeColor','none');
    hold on;
    
end

set(gca,'XTick',xbin,'XTickLabel',{'left','right','left','right'}, ...
    'YLim',[0 25],'LineWidth',1,'FontSize',18);

ylabel('number of significant sequences');
xlabel('(food restricted)         (water restricted)');
title(sprintf('%d rats, %d sessions',length(cfg.rats),length(fd)));

if cfg.writeFiles
     cd([fd{iFD},'\files']) %cd('files');
    base_fn = ExpKeys.goodSWR{1}(1:15);
    
    fn = cat(2,base_fn,'_allSessionsCorr.png');
    print(gcf,'-dpng','-r300',fn);
    
    cd ..
end
