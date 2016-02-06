% plot event triggered spectrograms



cfg.output_fd   = 'E:\Documents\TmazePaper\visuals'; % file directory for saving image files
cfg.output_fn   = 'EventTriggeredSpectrograms';
cfg.input_fd = 'E:\Documents\TmazePaper\data';
cfg.input_fn = 'TFRdata';
cfg.warningoff  = 1; % 1 don't show 8500 warnings that ft displays; else show warnings
cfg.plotOutput = 0;
cfg.plotRunning = 0; % if 1, plot the fake events too, if 0 don't

% must be the same as in ThreePSDs script
cfg.sessionFD   = {'E:\data\promoted\R042\R042-2013-08-16','E:\data\promoted\R044\R044-2013-12-20','E:\data\promoted\R050\R050-2014-04-03','E:\data\promoted\R064\R064-2015-04-23'};

cfg.rats = {'R042','R044','R050','R064'};

%%

iWasHere = pwd;
cd(cfg.input_fd);
load(cfg.input_fn)
position = 1;
m = length(rats);

if cfg.plotRunning;
    n = 2;
else
    n = 1;
end

for iRat = 1:length(cfg.rats)
    
    cd(cfg.sessionFD{iRat})
    LoadCandidates % for number of events
    
    [~,sessionID,~] = fileparts(pwd);
    
    
    TFR = TFRdata.(rats{iRat}).TFR;
    TFR_running = TFRdata.(rats{iRat}).TFR_running;
    
    % generate plot for candidates
    
    subplot(m,n,position)
   
    cfg_temp = [];
    cfg_temp.colorbar = 'no';
    
    cfg_temp.channel = TFR.label; 
    ft_singleplotTFR(cfg_temp, TFR);
    title([sessionID,': ',num2str(length(evt.tstart)),' candidates'], 'FontSize',8);
    grid on
    pbaspect(gca,[1 1 1])
    position = position + 1;
    
    % generate plot for running
    
    if cfg.plotRunning
        subplot(m,n,position)
        
        cfg_temp = []; 
        
        cfg_temp.channel = TFR_running.label; cfg_temp.colorbar = 'no';
        ft_singleplotTFR(cfg_temp, TFR_running);
        title([sessionID,': ',num2str(length(evt.tstart)),' running intervals'], 'FontSize',8);
        grid on
        pbaspect(gca,[1 1 1])
        position = position + 1;
    end
    
    
end

cd(cfg.output_fd)
% print(gcf,'-dpng','-r300',[cfg.output_fn,'.png']);
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);
print(gcf,'-dpdf','-r300',cfg.output_fn);
% print(gcf,'-depsc','-r300',cfg.output_fn);


cd(iWasHere)