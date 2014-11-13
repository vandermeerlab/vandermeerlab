function PlotTSDfromIV(cfg_in,iv,tsd_in)
% function PlotTSDfromIV(cfg,iv,tsd)
%
% create interval data from tsd by tresholding
%
% INPUTS:
%
% tsd_in: input tsd
%
% CFG OPTIONS:
% cfg.display = 'tsd'; % 'iv'
% cfg.mode = 'edges'; % 'center'
% cfg.width = 0.5; % in s
% cfg.subplotdim = [10 8];
% cfg.target = []; % if tsd has multiple data dimensions
%
%
% MvdM 2014-06-24

cfg_def.display = 'tsd'; % 'iv'
cfg_def.mode = 'edges'; % 'center'
cfg_def.width = 0.2; % in s
cfg_def.subplotdim = [10 8];
cfg_def.bgcol = 'k';
cfg_def.fgcol = 'r';

cfg = ProcessConfig2(cfg_def,cfg_in); % should take whatever is in cfg_in and put it into cfg!
mfun = mfilename;

% check if conditions are in place
nData = size(tsd_in.data,1);
if nData > 1
    if ~isempty(cfg.target)
        temp_data = getd(tsd_in,cfg.target);
    else
       error('Multiple data dimensions exist but no label is specified.'); 
    end
else
    temp_data = tsd_in.data;
end

% find indices for iv
switch cfg.mode
    case 'edges'
        tstart_idx = nearest_idx3(iv.tstart,tsd_in.tvec);
        tend_idx = nearest_idx3(iv.tend,tsd_in.tvec);
    case 'center'
        ctr = mean(cat(2,iv.tstart,iv.tend),2);
        tstart_idx = nearest_idx3(ctr-cfg.width/2,tsd_in.tvec);
        tend_idx = nearest_idx3(ctr+cfg.width/2,tsd_in.tvec);
end

switch cfg.display
    case 'tsd' % plot tsd with highlighted iv
        
        plot(tsd_in.tvec,temp_data,cfg.bgcol,'MarkerSize',1);
        hold on;
        
        for iI = 1:length(tstart_idx)
        
            plot(tsd_in.tvec(tstart_idx(iI):tend_idx(iI)),temp_data(tstart_idx(iI):tend_idx(iI)),cfg.fgcol,'MarkerSize',1);
            
        end
            
    case 'iv'
        
        nppf = prod(cfg.subplotdim);
        for iI = 1:length(tstart_idx)
            
            figno = ceil(iI./nppf);
            plotno = mod(iI-1,nppf) + 1; 
            maximize;
            
            figure(figno);
            subtightplot(cfg.subplotdim(1),cfg.subplotdim(2),plotno);
            
            plot(tsd_in.tvec(tstart_idx(iI):tend_idx(iI)),temp_data(tstart_idx(iI):tend_idx(iI)),cfg.fgcol,'MarkerSize',1);
            hold on;
            
            axis off; axis tight;
            
        end
    
end % switch cfg.display
