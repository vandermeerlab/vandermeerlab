function AnotherTCplot(cfg_in,S,TC,pos)
% function AnotherTCplot(cfg_in,S,TC,pos)
%
% plot some linearized tuning tuning curves with their scatterfields
%
% configs:
%
% cfg_def.npf = 6; % number of subplots per figure
% cfg_def.write_output = 0;

cfg_def.npf = 6; % number of subplots per figure
cfg_def.write_output = 0;

cfg = ProcessConfig2(cfg_def,cfg_in);

% should also create a version that plots left and right versions side by
% side in a 4-column, 2-row layout

for iC = 1:length(S.t)
   
    % find correct place to plot
    fno = floor((iC-1)/cfg.npf) + 1;
    sno = mod(iC,cfg.npf); if sno == 0, sno = cfg.npf; end
    
    fh = figure(fno); set(fh,'Color',[0 0 0]);
    
    % plot tuning curve
    subplot(2,6,1+(sno-1)*2);
    
    smoo_pfR{iC} = TC.tc(:,iC);
    
    ah = area(smoo_pfR{iC}); set(ah,'FaceColor',[0 0 0.5],'EdgeColor',[0 0 0.5]);
    axis off; set(gca,'YLim',[0 max(30,max(smoo_pfR{iC}))]); hold on;
    th = text(1,max(smoo_pfR{iC}),num2str(round(max(smoo_pfR{iC})))); set(th,'Color',[1 1 1]);
    
    % see if any fields were found, if so plot peaks also
    if ~isempty(TC.peak_loc{iC})
        pks.loc = TC.peak_loc{iC};
        pks.val = smoo_pfR{iC}(pks.loc);
        
        plot(pks.loc,pks.val,'.','MarkerSize',10,'Color',[0 0.7 0]);
    else
        set(ah,'FaceColor',[0.3 0.3 0.3],'EdgeColor',[0.3 0.3 0.3]);
    end
      
    % plot scatterfield
    subplot(2,6,sno*2);
    
    plot(getd(pos,'y'),getd(pos,'x'),'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
    hold on; axis off

    spk = S.t{iC};
    spk_x = interp1(pos.tvec,getd(pos,'x'),spk,'nearest');
    spk_y = interp1(pos.tvec,getd(pos,'y'),spk,'nearest');
    
    h = plot(spk_y,spk_x,'.r','MarkerSize',10);
    
    [~,fn,fe] = fileparts(S.label{iC});
    fs = cat(2,fn,fe);
    
    th = title(sprintf('c%d %s',iC,fs)); set(th,'Interpreter','none','Color',[1 1 1],'FontSize',8);
    maximize;
    
    if cfg.write_output & sno == cfg.npf
       fn_base = cat(2,S.cfg.SessionID,'-fields');
       fn = cat(2,fn_base,num2str(fno),'.png');
       set(gcf, 'InvertHardCopy', 'off');
       print(gcf,'-r300','-dpng',fn);
    end
    
end