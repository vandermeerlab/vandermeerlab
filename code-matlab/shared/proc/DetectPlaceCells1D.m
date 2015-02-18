function [idx,peak_idx,peak_loc] = DetectPlaceCells1D(cfg_in,tc)
% function [idx,peak_idx,peak_loc] = DetectPlaceCells1D(cfg,tc)
%
% inputs:
%
% tc: nBins x nCells matrix of tuning curves
%
% outputs:
%
% idx: indices (into tc) of detected place cells
% peak_idx: location of firing rate peak
% peak_loc: cell array of detected peak idx's
%
% cfg options:
%
%

cfg_def = [];
cfg_def.debug = 0;
cfg_def.p_thr = 7; % threshold for detecting place field, in Hz
cfg_def.p_thr_sd = 1; % in SD
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

for iC = 1:nCells

    this_tc = tc(:,iC);
    
    % first, decide if we are considering this cell at all based on average
    % firing rate -- could improve by having knowledge of whole-session
    % firing rate
    this_frate = nanmean(this_tc);

    if ~(max(this_tc) > cfg.p_thr & this_frate < cfg.max_meanfr) % too low or too high firing, reject
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
    
    %[pks.loc,pks.full] = find_fields(pf_z,cfg.p_thr_sd,'min_size',cfg.minSize);
    [pks.loc,pks.full] = find_fields(this_tc,cfg.p_thr,'min_size',cfg.minSize);
    
    if ~isempty(pks.loc) % some peaks were detected
        pks.val = this_tc(pks.loc); % get firing rate at peak
        
        keep = pks.val > cfg.p_thr;
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
        idx(ctr) = iC;
        [~,peak_idx(ctr)] = max(this_tc); % to determine where this cell appears in overall ordering of TC
        peak_loc{ctr} = pks.loc;
        
        ctr = ctr + 1;
    end % of detection
    
end

[~,sort_idx] = sort(peak_idx,'ascend');
idx = idx(sort_idx);
peak_idx = peak_idx(sort_idx);
peak_loc = peak_loc(sort_idx);
function [fields_out,fields_full] = find_fields(tc,thr,varargin)

min_size = 3;
extract_varargin;

% initialize
stack = find(tc >= thr);
wrap_point = length(tc);
fields_out = {}; fields_full = {}; fields = {};

% find fields
fieldCount = 0;
while ~isempty(stack)
    
    fieldCount = fieldCount + 1;
    item = stack(1);
    if length(stack) == 1
        stack = [];
    else
        stack = stack(2:end);
    end
    [stack,fields{fieldCount}] = growField(item,stack,wrap_point);
    
end

% process fields (output index of max firing rate for each field of sufficient size)
if fieldCount
    goodFieldCount = 1;
    for iF = 1:length(fields)
        if length(fields{iF}) >= min_size % size passed
            [val,ind] = max(tc(fields{iF}));
            fields_full{goodFieldCount} = fields{iF};
            fields_out{goodFieldCount} = fields{iF}(ind);
            goodFieldCount = goodFieldCount + 1;
        end
    end
end

fields_out = cell2mat(fields_out);

function [overall_stack,field_out] = growField(item_in,overall_stack,wrap_point)

step = 2;
wraparound = 0;
temp_stack = item_in; % start with input in the stack
field_items = []; % track places covered

while ~isempty(temp_stack);

    item = temp_stack(1);
    field_items = cat(1,field_items,item);
    if length(temp_stack) == 1
        temp_stack = [];
    else
        temp_stack = temp_stack(2:end);
    end
    
    % find sid_s close enough to current
    min_sid = item - step;
    if min_sid < 1
        min_sid = wrap_point + min_sid;
    end

    max_sid = item + step;
    if max_sid > wrap_point
        max_sid = max_sid - wrap_point;
    end

    if min_sid < max_sid

        temp_ind = overall_stack >= min_sid & overall_stack <= max_sid;

    else % wraparound

       wraparound = 1;
       temp_ind = overall_stack >= min_sid | overall_stack <= max_sid;
        
    end

    temp_stack = cat(1,temp_stack,overall_stack(temp_ind));
    overall_stack = overall_stack(~temp_ind);
    
end

field_out = field_items;