function [fields_out,fields_full] = find_fields(cfg_in,tc)
%% find_fields Find place fields
%   [fields_out,fields_full] = find_fields(cfg_in,tc) returns the peak indicies and
%   containing indicies of local maxima in a 1-D array. (used in DetectPlaceCells1D.m)

cfg_def = [];
cfg_def.thr = 5; %Threshold for firing rate of place cells in Hz
cfg_def.min_size = 3; %Minimum size of place field (units?)

cfg = ProcessConfig2(cfg_def,cfg_in);

% initialize
stack = find(tc >= cfg.thr);
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
        if length(fields{iF}) >=cfg.min_size % size passed
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