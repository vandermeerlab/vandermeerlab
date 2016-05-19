function h = PlotTSDfromIV(cfg_in,iv_in,tsd_in)
% function PlotTSDfromIV(cfg,iv_in,tsd_in)
%
% display intervals defined relative to tsd
%
% INPUTS:
%
% iv_in: input iv
% tsd_in: input tsd
%
% OUTPUTS:
%
% h: struct with handles to LFP and highlighted iv's
%
% CFG OPTIONS:
% cfg_def.display = 'tsd'; % {'tsd','iv'}, if 'tsd' then plot tsd and
%   highlight iv's within it; if 'iv' then only plot iv's
% cfg_def.bgcol = 'k'; % tsd color outside iv's
% cfg_def.fgcol = 'r'; % tsd color within iv's
% cfg_def.target = []; % which tsd channel to use if tsd has multiple data dimensions
%
% ONLY USED IN 'IV' DISPLAY MODE:
% cfg_def.mode = 'edges'; % {'edges','center'}: plot each iv either using a
%   specific time window ('center') or relative to each iv's edges ('edges')
% cfg_def.width = 0.2; % size of time window in s (for 'center') or time
%   window to be added to each edge (for 'edges')
% cfg_def.subplotdim = [10 8]; % for 'iv' display mode only, specifies
%   subplots to use in single figure
% cfg_def.title = []; % for each iv, add contents of specified usr
%  field as title (SHOULD GENERALIZE TO ARBITRARY STRING WITH MULTIPLE
%  FIELDS, FORMATTING...)
%
% MvdM 2014-06-24, edit 2016-01-0y to add title option

cfg_def = [];
cfg_def.verbose = 1;
cfg_def.display = 'tsd'; % 'iv'
cfg_def.mode = 'edges'; % 'center'
cfg_def.width = 0.2; % in s
cfg_def.subplotdim = [10 8];
cfg_def.bgcol = 'k';
cfg_def.fgcol = 'r';
cfg_def.MarkerSize = 10;
cfg_def.iv_only = 0; % if 1, don't plot tsd
cfg_def.title = [];

cfg = ProcessConfig(cfg_def,cfg_in); % should take whatever is in cfg_in and put it into cfg!
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
tstart_idx = nearest_idx3(iv_in.tstart,tsd_in.tvec);
tend_idx = nearest_idx3(iv_in.tend,tsd_in.tvec);

% if iv display mode, also find indices for background
switch cfg.mode
    case 'edges'
        bg_tstart_idx = nearest_idx3(iv_in.tstart-cfg.width,tsd_in.tvec);
        bg_tend_idx = nearest_idx3(iv_in.tend+cfg.width,tsd_in.tvec);
    case 'center'
        ctr = mean(cat(2,iv_in.tstart,iv_in.tend),2);
        bg_tstart_idx = nearest_idx3(ctr-cfg.width/2,tsd_in.tvec);
        bg_tend_idx = nearest_idx3(ctr+cfg.width/2,tsd_in.tvec);
end

switch cfg.display
    case 'tsd' % plot tsd with highlighted iv
        
        if ~cfg.iv_only
        h.LFP = plot(tsd_in.tvec,temp_data,cfg.bgcol,'MarkerSize',cfg.MarkerSize);
        end
        hold on;
        
        h.LFP_iv = nan(size(tstart_idx));
        for iI = 1:length(tstart_idx)
        
            h.LFP_iv(iI) = plot(tsd_in.tvec(tstart_idx(iI):tend_idx(iI)),temp_data(tstart_idx(iI):tend_idx(iI)),cfg.fgcol,'MarkerSize',cfg.MarkerSize);
            
        end
            
    case 'iv'
        
        nppf = prod(cfg.subplotdim);
        for iI = 1:length(tstart_idx)
            
            figno = ceil(iI./nppf);
            plotno = mod(iI-1,nppf) + 1; 
            maximize;
            
            figure(figno);
            subtightplot(cfg.subplotdim(1),cfg.subplotdim(2),plotno);
            
            h.LFP(iI) = plot(tsd_in.tvec(bg_tstart_idx(iI):bg_tend_idx(iI)),temp_data(bg_tstart_idx(iI):bg_tend_idx(iI)),cfg.bgcol,'MarkerSize',1);
            hold on;
            h.LFP_iv(iI) = plot(tsd_in.tvec(tstart_idx(iI):tend_idx(iI)),temp_data(tstart_idx(iI):tend_idx(iI)),cfg.fgcol,'MarkerSize',1);
                        
            axis off; axis tight;
            if ~isempty(cfg.title)
               all_usr = iv_in.usr.(cfg.title);
               title(all_usr(iI));
            end
            
        end
    
end % switch cfg.display
