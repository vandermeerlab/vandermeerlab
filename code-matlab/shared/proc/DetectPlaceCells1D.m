function tc_out = DetectPlaceCells1D(cfg_in,tc)
% function [idx,peak_idx,peak_loc] = DetectPlaceCells1D(cfg,tc)
%
% inputs:
%
% tc: nBins x nCells matrix of tuning curves
%
% outputs:
%
% tc_out.tc: original nCells x nBins tuning curves
% tc_out.template_idx: indices (into tc) of detected place cells
% tc_out.peak_idx: location of firing rate peak
% tc_out.peak_loc: cell array of detected peak idx's
% tc_out.field_loc = fpeak_val; single field locations (in cm?)
% tc_out.peak_idx = peak_idx; multiple field indicies
% tc_out.peak_loc = peak_loc; multiple field locations
%
%
% cfg options:
%
%

cfg_def = [];
cfg_def.debug = 0;
cfg_def.thr = 7; % threshold for detecting place field, in Hz
cfg_def.thr_sd = 1; % in SD
cfg_def.mean_norm = 0.33; % max mean firing rate (peak is 1)
cfg_def.minSize = 4;
%cfg_def.maxSize = 20;
cfg_def.max_meanfr = 10;
cfg_def.S = []; % optionally provide actual spike data for computing true mean frate, nSpikes per field etc.
cfg_def.pos = []; % optionally provide position data for getting nSpikes per field
cfg_def.nSpikesInField = [];

cfg = ProcessConfig2(cfg_def,cfg_in);

nCells = size(tc,2);

ctr = 1; % counter for detected cells

% !!!!!!!!!!
tc = tc'; % DO A BETTER FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!

peak_idx = []; peak_loc = []; template_idx = [];
for iC = 1:nCells

    this_tc = tc(:,iC);
    
    % first, decide if we are considering this cell at all based on average
    % firing rate -- could improve by having knowledge of whole-session
    % firing rate
    this_frate = nanmean(this_tc);

    if ~(max(this_tc) > cfg.thr & this_frate < cfg.max_meanfr) % too low or too high firing, reject
        if cfg.debug, fprintf('Cell %d rejected: max TC %.1f, mean FR %.1f\n',iC,max(this_tc),this_frate); end
        continue;
    end
    
    pf_z = zscore(this_tc);
    pf_norm = this_tc./max(this_tc);
    mean_norm = mean(pf_norm);
    %
    if mean_norm > cfg.mean_norm
        if cfg.debug, fprintf('Cell %d rejected: mean norm %.2f\n',iC,mean_norm); end
        continue;
    end
    
    %[pks.loc,pks.full] = find_fields(pf_z,cfg.thr_sd,'min_size',cfg.minSize);
%     [pks.loc,pks.full] = find_fields(this_tc,cfg.thr,'min_size',cfg.minSize);
    [pks.loc,pks.full] = find_fields(cfg,this_tc);
        
    if ~isempty(pks.loc) % some peaks were detected
        pks.val = this_tc(pks.loc); % get firing rate at peak
        
        keep = pks.val > cfg.thr;
        pks.loc = pks.loc(keep); pks.val = pks.val(keep); pks.full = pks.full(keep);
        if cfg.debug & isempty(keep), fprintf('Cell %d, no peaks passed threshold.\n',iC); end
    else
        if cfg.debug, fprintf('Cell %d, no peaks detected.\n',iC); end
    end
    
    if ~isempty(pks.loc) & ~isempty(cfg.nSpikesInField) % get number of spikes in field
        
        clear nSpikes;
        for iF = 1:length(pks.loc)
            
           field_start = pks.full{iF}(1)-1; % kind of a hack...
           field_end = pks.full{iF}(end)+1; % could instead requre peak > x, edges at y
           
           % get field entries and exits as an iv
           cfg_infield = [];
           cfg_infield.method = 'raw'; cfg_infield.threshold = [field_start field_end]; 
           cfg_infield.dcn = 'range'; cfg_infield.target = 'z';
           infield_iv = TSDtoIV(cfg_infield,cfg.pos);
           
           % now get number of spikes
           this_s = restrict(cfg.S,infield_iv);
           nSpikes(iF) = length(this_s.t{iC});
            
           if cfg.debug & nSpikes(iF) < cfg.nSpikesInField, fprintf('Cell %d, field %d rejected (nSpikes = %d).\n',iC,iF,nSpikes(iF)); end
           
        end
        
        keep = nSpikes > cfg.nSpikesInField;
        pks.loc = pks.loc(keep); pks.val = pks.val(keep); pks.full = pks.full(keep);
        
    end
        
    if ~isempty(pks.loc) %  some fields remain, process
        template_idx(ctr) = iC;
        [~,peak_idx(ctr)] = max(this_tc); % to determine where this cell appears in overall ordering of TC
        peak_loc{ctr} = pks.loc;
        
        ctr = ctr + 1;
    end % of detection
    
end

[~,sort_idx] = sort(peak_idx,'ascend');
template_idx = template_idx(sort_idx);
peak_idx = peak_idx(sort_idx);
peak_loc = peak_loc(sort_idx);

% create field template
fpeak_val = []; fpeak_idx = [];
iPeak = 1;
for iP = 1:length(peak_loc)
    fpeak_val(iPeak:iPeak+length(peak_loc{iP})-1) = peak_loc{iP};
    fpeak_idx(iPeak:iPeak+length(peak_loc{iP})-1) = template_idx(iP);
    
    iPeak = iPeak + length(peak_loc{iP});
end
[fpeak_val,sort_idx] = sort(fpeak_val,'ascend');
fpeak_idx = fpeak_idx(sort_idx);

%store all data in a single struct
tc_out.tc = tc;
tc_out.template_idx = template_idx;
tc_out.field_template_idx = fpeak_idx;
tc_out.field_loc = fpeak_val;
tc_out.peak_idx = peak_idx;
tc_out.peak_loc = peak_loc;


