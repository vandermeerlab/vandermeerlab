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
cfg_def.p_thr = 5; % threshold for detecting place field, in Hz
cfg_def.p_thr_sd = 2; % in SD
cfg_def.max_meanfr = 2;
cfg_def.S = []; % optionally provide actual spike data for computing true mean frate, nSpikes per field etc.
cfg_def.pos = []; % optionally provide position data for getting nSpikes per field
cfg_def.nSpikesInField = [];

cfg = ProcessConfig2(cfg_def,cfg_in);

nCells = size(tc,2);

ctr = 1; % counter for detected cells

for iC = 1:nCells

    this_tc = tc(:,iC);
    
    % first, decide if we are considering this cell at all based on average
    % firing rate
    if ~isempty(cfg.S) % we have spike data to get good estimate
        this_frate = length(S.t{iC})./(S.cfg.ExpKeys.postrecord(end)-S.cfg.ExpKeys.prerecord(1));
    else
        this_frate = nanmean(this_tc);
    end
    
    if ~(max(this_tc) > cfg.p_thr & this_frate < cfg.max_meanfr) % too low or too high firing, reject
        continue;
    end
    
    pf_z = (this_tc-nanmean(this_tc))./nanstd(this_tc);
    
    [pks.loc,pks.full] = find_fields(pf_z,cfg.p_thr_sd,'min_size',15);
    
    if ~isempty(pks.loc) % some peaks were detected
        pks.val = this_tc(pks.loc); % get firing rate at peak
        
        keep = pks.val > cfg.p_thr;
        pks.loc = pks.loc(keep); pks.val = pks.val(keep); pks.full = pks.full(keep);
    end
    
    if ~isempty(pks.loc) & ~isempty(cfg.nSpikesInField) % get number of spikes in field
        
        for iF = 1:length(pks.loc)
            
           field_start = pks.full{iF}(1);
           field_end = pks.full{iF}(end);
           
           % get field entries and exits as an iv
           cfg_infield = [];
           cfg_infield.method = 'raw'; cfg_infield.threshold = [field_start field_end]; 
           cfg_infield.dcn = 'range'; cfg_infield.target = 'z';
           infield_iv = TSDtoIV(cfg_infield,cfg.pos);
           
           % now get number of spikes
            
        end
        
    end
        
    if ~isempty(pks.loc) % field still lives, process
        idx(ctr) = iC;
        [~,peak_idx(ctr)] = max(this_tc);
        peak_loc{ctr} = pks.loc;
        
        ctr = ctr + 1;
    end % of detection
    
end

[~,sort_idx] = sort(peak_idx,'ascend');
idx = idx(sort_idx);
peak_idx = peak_idx(sort_idx);
peak_loc = peak_loc(sort_idx);

function fields = find_fields(tc,thr,varargin)

min_size = 3;
extract_varargin;

% initialize
stack = find(tc >= thr);
wrap_point = length(tc);
fields = {};

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
        if length(fields{iF} >= min_size)
            [val,ind] = max(tc(fields{iF}));
            fields{goodFieldCount} = fields{iF}(ind);
            goodFieldCount = goodFieldCount + 1;
        end
    end
end

fields = cell2mat(fields);

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