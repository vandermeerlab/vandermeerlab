function tc_out = DetectPlaceCells1D(cfg_in,tc_in)
% function [idx,peak_idx,peak_loc] = DetectPlaceCells1D(cfg,tc)
%
% inputs:
%
% tc: tuning curves
%   - nBins x nCells matrix of tuning curves OR tc struct with tc field
%
% outputs:
%
% tc_out.tc: original nCells x nBins tuning curves
% tc_out.template_idx: 1 x nPlaceCells idxs (into tc) of detected place
%   cells, ordered by peak_idx
% tc_out.peak_idx: 1 x nPlaceCells idxs (into tc) indicating location of firing rate peak
% tc_out.peak_loc: 1 x nPlaceCells cell array of detected peak idxs (because one field can
%   have multiple peaks)
% tc_out.field_loc: 1 x nFields location of field peaks
% tc_out.peak_idx: ??
%
% cfg options:
% cfg_def.debug = 0; % set to 1 to get verbose output
% cfg_def.thr = 5; % min threshold for detecting candidate place fields in Hz
% cfg_def.minSize = 4; % minimum size (in bins) that must exceed threshold
% cfg_def.max_meanfr = 10; % max mean firing rate (exclude interneurons)
% cfg_def.nSpikesInField = []; % minimum number of spikes in field; if not
% empty, must also provide cfg_in.S and cfg.pos
% cfg_def.S = [];
% cfg_def.pos = [];
%
% MvdM, youkitan

if isstruct(tc_in)
    tc = tc_in.tc';
else
    tc = tc_in'; % FIX!!!!!
end

cfg_def = [];
cfg_def.debug = 0; % set to 1 to get verbose output
cfg_def.thr = 5; % min threshold for detecting candidate place fields in Hz
cfg_def.minSize = 4; % minimum size (in bins) that must exceed threshold
cfg_def.max_meanfr = 10; % max mean firing rate (exclude interneurons)
cfg_def.nSpikesInField = []; % minimum number of spikes in field; if not empty, must also provide cfg_in.S and cfg_in.pos
cfg_def.S = [];
cfg_def.pos = [];
cfg_def.verbose = 1; % controls whether ProcessConfig will warn you if there's a typo
mfun = mfilename;

cfg = ProcessConfig(cfg_def,cfg_in,mfun);

nCells = size(tc,2);

ctr = 1; % counter for detected cells

peak_idx = []; peak_loc = []; template_idx = [];

% get firing rate
if isstruct(tc_in) && isfield(tc_in,'usr') && isfield(tc_in.usr,'rate')
    % whole session firing rate
    frate = tc_in.usr.rate;
else
    frate = nanmean(tc);
end

for iC = 1:nCells

    this_tc = tc(:,iC);
    
    % first, decide if we are considering this cell at all based on average
    % firing rate
    
    this_frate = frate(iC);
    
    if ~(max(this_tc) > cfg.thr && this_frate < cfg.max_meanfr) % too low or too high firing, reject
        if cfg.debug, fprintf('Cell %d rejected: max TC %.1f, mean FR %.1f\n',iC,max(this_tc),this_frate); end
        continue;
    end
    
    cfg_ff = [];
    cfg_ff.thr = cfg.thr;
    cfg_ff.minSize = cfg.minSize;
    [pks.loc,pks.full] = find_fields(cfg_ff,this_tc);
    
    if ~isempty(pks.loc) % some peaks were detected
        pks.val = this_tc(pks.loc); % get firing rate at peak
        
        keep = pks.val > cfg.thr;
        pks.loc = pks.loc(keep); pks.val = pks.val(keep); pks.full = pks.full(keep);
        if cfg.debug && isempty(keep), fprintf('Cell %d, no peaks passed threshold.\n',iC); end
    else
        if cfg.debug, fprintf('Cell %d, no peaks detected.\n',iC); end
    end
    
    if ~isempty(pks.loc) && ~isempty(cfg.nSpikesInField) % get number of spikes in field
        
        clear nSpikes;
        for iF = 1:length(pks.loc)
            
           field_start = pks.full{iF}(1)-1; % kind of a hack...
           field_end = pks.full{iF}(end)+1; % could instead requre peak > x, edges at y
           
           % get field entries and exits as an iv
           cfg_infield = [];
           cfg_infield.method = 'raw'; cfg_infield.threshold = [field_start field_end]; 
           cfg_infield.dcn = 'range'; cfg_infield.target = 'z'; cfg_infield.verbose = 0;
           infield_iv = TSDtoIV(cfg_infield,cfg.pos);
           
           % now get number of spikes
           this_s = restrict(cfg.S,infield_iv);
           nSpikes(iF) = length(this_s.t{iC});
            
           if cfg.debug && nSpikes(iF) < cfg.nSpikesInField, fprintf('Cell %d, field %d rejected (nSpikes = %d).\n',iC,iF,nSpikes(iF)); end
           
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

% write history
if isstruct(tc_in) && isfield(tc_in,'cfg')
    tc_out.cfg.history = tc_in.cfg.history;
end
tc_out = History(tc_out,mfun,cfg);

end

