function PlotFields(cfg_in,tc,fields,pos,S)
% function PlotFields(cfg_in,tc,fields,pos,S)
%
% plots output of FindFields() superimposed on tuning curves alongside
% scatterfields, for inspection of detected fields
%
% INPUTS:
% cfg_def.npf = 6; % number of subplots per figure
% cfg_def.subsample = 10; % subsampling factor for position data
%
% fields: output of FindFields()
% pos: position data
% S: spikes (note idxs of fields should line up with this)
%
% OUTPUTS:
% (none)
%
% MvdM 2015-10-01 initial version

cfg_def.npf = 6;
cfg_def.subsample = 10;

cfg = ProcessConfig2(cfg_def,cfg_in);

for iC = 1:length(S.t)
    
    % find current figure number
    fno = floor((iC-1)/cfg.npf) + 1;
    sno = mod(iC,cfg.npf); if sno == 0, sno = cfg.npf; end
    
    fh = figure(fno); set(fh,'Color',[0 0 0]);
    
    % plot tuning curve
    subplot(2,6,1+(sno-1)*2);
    
    ah = area(tc(iC,:)); set(ah,'FaceColor',[0 0 0.7],'EdgeColor',[0 0 0.7]);
    axis off; set(gca,'YLim',[0 max(30,max(tc(iC,:)))]); hold on;
    th = text(1,max(tc(iC,:)),num2str(round(max(tc(iC,:))))); set(th,'Color',[1 1 1]);
    
    % plot fields if detected
    field_idx = find(fields.template_idx == iC);
    if isempty(field_idx)
        set(ah,'FaceColor',[0.3 0.3 0.3],'EdgeColor',[0.3 0.3 0.3]);
    else
        these_fields = fields.peak_loc{field_idx};
        plot(these_fields,tc(iC,these_fields),'.','MarkerSize',10,'Color',[0 0.7 0]);
    end
    
    % plot scatterfield
    subplot(2,6,sno*2);
    
    x = getd(pos,'x'); x = x(1:cfg.subsample:end);
    y = getd(pos,'y'); y = y(1:cfg.subsample:end);
    
    plot(y,x,'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
    hold on; axis off

    spk = S.t{iC};
    spk_x = interp1(pos.tvec,getd(pos,'x'),spk,'nearest');
    spk_y = interp1(pos.tvec,getd(pos,'y'),spk,'nearest');
    
    h = plot(spk_y,spk_x,'.r','MarkerSize',10);
    
    [fp fn fe] = fileparts(S.label{iC});
    fs = cat(2,fn,fe);
    
    th = title(sprintf('c%d %s',iC,fs)); set(th,'Interpreter','none','Color',[1 1 1],'FontSize',8);
    maximize;
      
end % of loop over cells