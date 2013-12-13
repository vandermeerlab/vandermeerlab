function [trl,event] = ft_trialfun_lineartracktone2(cfg,varargin)
% function [trl,event] = ft_trialfun_lineartracktone2(cfg,varargin)
%
%cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
%cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
%cfg.trialdef.block = 'both'; % could be 'value', 'risk'
%cfg.trialdef.cue = {'c1'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'}
%
% varargins:
% OutputFormat = 'Timestamps'; % could be 'Time'
% LeftRightMode = 'Position'; % could be 'PrevNosepoke' (method to determine L/R trial)
% EvtOnly = 0; % 
%
% MvdM 2013-07-03

OutputFormat = 'Timestamps'; % could be 'Time'
LeftRightMode = 'Position'; % could be 'PrevNosepoke' (method to determine L/R trial)
EvtOnly = 0; % don't segment into trials
debug_mode = 0;

extract_varargin;

%%
fn = FindFile('*Events.nev');

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

%% determine event strings to detect

% first check if an unused photobeam is on (from Rob's task, which ran in
% parallel with Sushant's LinearTrackTone)

idx = strmatch('TTL Input on AcqSystem1_0 board 0 port 1 value (0x0004)',EventStrings);

if ~isempty(idx)
    fprintf('NOTICE: photobeam 3 detected. Verifying data consistency...\n');
    
    idx2 = strmatch('TTL Input on AcqSystem1_0 board 0 port 1 value (0x0001)',EventStrings);
    idx3 = strmatch('TTL Input on AcqSystem1_0 board 0 port 1 value (0x0002)',EventStrings);
    
    if ~isempty(idx2) | ~isempty(idx3)
        error('0x0001 or 0x0002 events found; should not occur if photobeam 3 is on!');
    else
        fprintf('NOTICE: events OK, no 0x0001 or 0x0002 events found.\n');
    end
    
    pboff_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0004)';
    pb0_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0005)';
    pb1_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0006)';
else
    pboff_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000)';
    pb0_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0001)';
    pb1_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0002)';
end

foff_string = 'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0000)';
f0_string = 'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0001)';
f1_string = 'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0002)';

c1_string = '1 pellet cue';
c3_string = '3 pellet cue';
c5_string = '5 pellet cue';
lo_string = '2 or 4 pellet cue';
hi_string = '1 or 5 pellet cue';

%% NOTE ALL PROCESSING HAPPENS IN INDEX SPACE (INTO EVENTSTRINGS, MATCHING EVTIMESTAMPS)

%% step 1: collect all raw event onsets
pb0_idx = strmatch(pb0_string,EventStrings);
pb1_idx = strmatch(pb1_string,EventStrings);

f0_idx = strmatch(f0_string,EventStrings);
f1_idx = strmatch(f1_string,EventStrings);

c1_idx = strmatch(c1_string,EventStrings);
c3_idx = strmatch(c3_string,EventStrings);
c5_idx = strmatch(c5_string,EventStrings);
lo_idx = strmatch(lo_string,EventStrings);
hi_idx = strmatch(hi_string,EventStrings);

all_cue_idx = cat(1,c1_idx,c3_idx,c5_idx,lo_idx,hi_idx);
%% step 2: find matching offsets
all_pboff_idx = strmatch(pboff_string,EventStrings);

for iE = length(pb0_idx):-1:1
    offset = all_pboff_idx - pb0_idx(iE);
    offset(offset < 0) = Inf;
    [~,following_idx] = min(offset);
    pb0_off_idx(iE) = all_pboff_idx(following_idx);
end

for iE = length(pb1_idx):-1:1
    offset = all_pboff_idx - pb1_idx(iE);
    offset(offset < 0) = Inf;
    [~,following_idx] = min(offset);
    pb1_off_idx(iE) = all_pboff_idx(following_idx);
end

%
all_foff_idx = strmatch(foff_string,EventStrings);

for iE = length(f0_idx):-1:1
    offset = all_foff_idx - f0_idx(iE);
    offset(offset < 0) = Inf;
    [~,following_idx] = min(offset);
    f0_off_idx(iE) = all_foff_idx(following_idx);
end

for iE = length(f1_idx):-1:1
    offset = all_foff_idx - f1_idx(iE);
    offset(offset < 0) = Inf;
    [~,following_idx] = min(offset);
    f1_off_idx(iE) = all_foff_idx(following_idx);
end


%% step 3: find correct nosepokes (feeder is next), preceding cue, and latencies
iCorrect = 1; clear pb0_correct_idx pb1_correct_idx cue0_correct_idx cue1_correct_idx;
for iE = 1:length(pb0_idx)
    
    curr_event_idx = pb0_idx(iE);
    next_event = EventStrings(curr_event_idx+1);
    
    if strcmp(next_event,'Feeder 0 nosepoke') | strcmp(next_event,f0_string)
        pb0_correct_idx(iCorrect) = curr_event_idx;
        
        % find preceding cue
        
        % first need preceding feeder fire
        preceding_f = cat(1,f0_idx,f1_idx);
        preceding_f = preceding_f(preceding_f < curr_event_idx);
        if ~isempty(preceding_f)
            preceding_f = max(preceding_f);
        else
            preceding_f = 1;
        end
        
        preceding_events_idx = preceding_f:curr_event_idx-1;
        preceding_cue_idx = intersect(preceding_events_idx,all_cue_idx);
        
        if ~isempty(preceding_cue_idx)
        
        cue0_correct_idx(iCorrect) = preceding_cue_idx(end);
        
        % latency
        cue0_correct_lat(iCorrect) = EVTimeStamps(cue0_correct_idx(iCorrect))-EVTimeStamps(pb0_correct_idx(iCorrect));
        cue0_correct_lat(iCorrect) = cue0_correct_lat(iCorrect)*10^-6;
        
        else % should not happen -- no preceding cue
            fprintf('WARNING: no preceding cue found for event %d\n',iE);
            cue0_correct_idx(iCorrect) = NaN;
            cue0_correct_lat(iCorrect) = NaN;           
        end
        
        iCorrect = iCorrect + 1;
    end
    
end
fprintf('f0: %d correct pokes found (%d feeder fires)\n',iCorrect-1,length(f0_idx));

% if number of pokes and fires is not equal, correct
if iCorrect-1 ~= length(f0_idx)
   fprintf('WARNING: inserting missing events\n');
   [pb0_correct_idx,f0_idx] = equalize_pbf(pb0_correct_idx,f0_idx);
end

iCorrect = 1;
for iE = 1:length(pb1_idx)
    
    curr_event_idx = pb1_idx(iE);
    next_event = EventStrings(curr_event_idx+1);
    
    if strcmp(next_event,'Feeder 1 nosepoke') | strcmp(next_event,f1_string)
        pb1_correct_idx(iCorrect) = curr_event_idx;
        
        % find preceding cue
        
        % first need preceding feeder fire
        preceding_f = cat(1,f0_idx,f1_idx);
        preceding_f = preceding_f(preceding_f < curr_event_idx);
        if ~isempty(preceding_f)
            preceding_f = max(preceding_f);
        else
            preceding_f = 1;
        end
        
        preceding_events_idx = preceding_f:curr_event_idx-1;
        preceding_cue_idx = intersect(preceding_events_idx,all_cue_idx);
        
        if ~isempty(preceding_cue_idx)
            
            cue1_correct_idx(iCorrect) = preceding_cue_idx(end);
            
            % latency
            cue1_correct_lat(iCorrect) = EVTimeStamps(cue1_correct_idx(iCorrect))-EVTimeStamps(pb1_correct_idx(iCorrect));
            cue1_correct_lat(iCorrect) = cue1_correct_lat(iCorrect)*10^-6;
            
        else % should not happen -- no preceding cue
            fprintf('WARNING: no preceding cue found for event %d\n',iE);
            cue0_correct_idx(iCorrect) = NaN;
            cue0_correct_lat(iCorrect) = NaN;
        end
        
        iCorrect = iCorrect + 1;
    end
    
end
fprintf('f1: %d correct pokes found (%d feeder fires)\n',iCorrect-1,length(f1_idx));

% if number of pokes and fires is not equal, correct
if iCorrect-1 ~= length(f1_idx)
   fprintf('WARNING: inserting missing events\n');
   [pb1_correct_idx,f1_idx] = equalize_pbf(pb1_correct_idx,f1_idx);
end

cue_correct_idx = cat(2,cue0_correct_idx,cue1_correct_idx);
cue_correct_lat = cat(2,cue0_correct_lat,cue1_correct_lat);

%% step 4: find left/right trial for cues (based on last nosepoke)

switch LeftRightMode
    case 'Nosepoke'
        
        % c1
        iLeft = 1; iRight = 1; clear c1_idx_left c1_idx_right
        for iE = 1:length(c1_idx)
            
            curr_event_idx = c1_idx(iE);
            preceding_events_idx = 1:curr_event_idx-1;
            
            preceding_pb0_idx = intersect(preceding_events_idx,pb0_idx);
            preceding_pb1_idx = intersect(preceding_events_idx,pb1_idx);
            
            if isempty(preceding_pb0_idx) & isempty(preceding_pb1_idx)
                error('No preceding nosepoke, cannot determine trial direction.');
            elseif isempty(preceding_pb0_idx)
                preceding_pb0_idx = 0;
            elseif isempty(preceding_pb1_idx)
                preceding_pb1_idx = 0;
            end
            
            if max(preceding_pb0_idx) < max(preceding_pb1_idx)
                c1_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                c1_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
            
        end
        
        % c3
        iLeft = 1; iRight = 1; clear c3_idx_left c3_idx_right
        for iE = 1:length(c3_idx)
            
            curr_event_idx = c3_idx(iE);
            preceding_events_idx = 1:curr_event_idx-1;
            
            preceding_pb0_idx = intersect(preceding_events_idx,pb0_idx);
            preceding_pb1_idx = intersect(preceding_events_idx,pb1_idx);
            
            if isempty(preceding_pb0_idx) & isempty(preceding_pb1_idx)
                error('No preceding nosepoke, cannot determine trial direction.');
            elseif isempty(preceding_pb0_idx)
                preceding_pb0_idx = 0;
            elseif isempty(preceding_pb1_idx)
                preceding_pb1_idx = 0;
            end
            
            if max(preceding_pb0_idx) < max(preceding_pb1_idx)
                c3_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                c3_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
            
        end
        
        % c5
        iLeft = 1; iRight = 1; clear c5_idx_left c5_idx_right
        for iE = 1:length(c5_idx)
            
            curr_event_idx = c5_idx(iE);
            preceding_events_idx = 1:curr_event_idx-1;
            
            preceding_pb0_idx = intersect(preceding_events_idx,pb0_idx);
            preceding_pb1_idx = intersect(preceding_events_idx,pb1_idx);
            
            if isempty(preceding_pb0_idx) & isempty(preceding_pb1_idx)
                error('No preceding nosepoke, cannot determine trial direction.');
            elseif isempty(preceding_pb0_idx)
                preceding_pb0_idx = 0;
            elseif isempty(preceding_pb1_idx)
                preceding_pb1_idx = 0;
            end
            
            if max(preceding_pb0_idx) < max(preceding_pb1_idx)
                c5_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                c5_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
            
        end
        
        % lo
        iLeft = 1; iRight = 1; clear lo_idx_left lo_idx_right
        for iE = 1:length(lo_idx)
            
            curr_event_idx = lo_idx(iE);
            preceding_events_idx = 1:curr_event_idx-1;
            
            preceding_pb0_idx = intersect(preceding_events_idx,pb0_idx);
            preceding_pb1_idx = intersect(preceding_events_idx,pb1_idx);
            
            if isempty(preceding_pb0_idx) & isempty(preceding_pb1_idx)
                error('No preceding nosepoke, cannot determine trial direction.');
            elseif isempty(preceding_pb0_idx)
                preceding_pb0_idx = 0;
            elseif isempty(preceding_pb1_idx)
                preceding_pb1_idx = 0;
            end
            
            if max(preceding_pb0_idx) < max(preceding_pb1_idx)
                lo_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                lo_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
            
        end
        
        % hi
        iLeft = 1; iRight = 1; clear hi_idx_left hi_idx_right
        for iE = 1:length(hi_idx)
            
            curr_event_idx = hi_idx(iE);
            preceding_events_idx = 1:curr_event_idx-1;
            
            preceding_pb0_idx = intersect(preceding_events_idx,pb0_idx);
            preceding_pb1_idx = intersect(preceding_events_idx,pb1_idx);
            
            if isempty(preceding_pb0_idx) & isempty(preceding_pb1_idx)
                error('No preceding nosepoke, cannot determine trial direction.');
            elseif isempty(preceding_pb0_idx)
                preceding_pb0_idx = 0;
            elseif isempty(preceding_pb1_idx)
                preceding_pb1_idx = 0;
            end
            
            if max(preceding_pb0_idx) < max(preceding_pb1_idx)
                hi_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                hi_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
            
        end
        
        all_cueleft_idx = cat(2,c1_idx_left,c3_idx_left,c5_idx_left,lo_idx_left,hi_idx_left);
        all_cueright_idx = cat(2,c1_idx_right,c3_idx_right,c5_idx_right,lo_idx_right,hi_idx_right);
        all_cue_idx = cat(2,all_cueleft_idx,all_cueright_idx);
               
    case 'Position'
        
        load(FindFile('*vt.mat'))
        X_REF = 320; % reference point (midway along track) to decide left/right
        
        iLeft = 1; iRight = 1; clear c1_idx_left c1_idx_right
        for iE = 1:length(c1_idx)
            
            curr_event_idx = c1_idx(iE);
            
            curr_ts = EVTimeStamps(curr_event_idx);
            curr_ts = curr_ts*10^-6;
            
            median_x = median(interp1(Range(x),Data(x),curr_ts-1:0.1:curr_ts));
            
            if median_x > X_REF 
                c1_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                c1_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
        end
        
        iLeft = 1; iRight = 1; clear c3_idx_left c3_idx_right
        for iE = 1:length(c3_idx)
            
            curr_event_idx = c3_idx(iE);
            
            curr_ts = EVTimeStamps(curr_event_idx);
            curr_ts = curr_ts*10^-6;
            
            median_x = median(interp1(Range(x),Data(x),curr_ts-1:0.1:curr_ts));
            
            if median_x > X_REF 
                c3_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                c3_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
        end
        
        iLeft = 1; iRight = 1; clear c5_idx_left c5_idx_right
        for iE = 1:length(c5_idx)
            
            curr_event_idx = c5_idx(iE);
            
            curr_ts = EVTimeStamps(curr_event_idx);
            curr_ts = curr_ts*10^-6;
            
            median_x = median(interp1(Range(x),Data(x),curr_ts-1:0.1:curr_ts));
            
            if median_x > X_REF 
                c5_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                c5_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
        end
        
        iLeft = 1; iRight = 1; clear lo_idx_left lo_idx_right
        for iE = 1:length(lo_idx)
            
            curr_event_idx = lo_idx(iE);
            
            curr_ts = EVTimeStamps(curr_event_idx);
            curr_ts = curr_ts*10^-6;
            
            median_x = median(interp1(Range(x),Data(x),curr_ts-1:0.1:curr_ts));
            
            if median_x > X_REF 
                lo_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                lo_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
        end
        
        iLeft = 1; iRight = 1; clear hi_idx_left hi_idx_right
        for iE = 1:length(hi_idx)
            
            curr_event_idx = hi_idx(iE);
            
            curr_ts = EVTimeStamps(curr_event_idx);
            curr_ts = curr_ts*10^-6;
            
            median_x = median(interp1(Range(x),Data(x),curr_ts-1:0.1:curr_ts));
            
            if median_x > X_REF 
                hi_idx_left(iLeft) = curr_event_idx;
                iLeft = iLeft + 1;
            else
                hi_idx_right(iRight) = curr_event_idx;
                iRight = iRight + 1;
            end
            
        end
        
        %%% add risk versions here too...
        
        all_cueleft_idx = cat(2,c1_idx_left,c3_idx_left,c5_idx_left,lo_idx_left,hi_idx_left);
        all_cueright_idx = cat(2,c1_idx_right,c3_idx_right,c5_idx_right,lo_idx_right,hi_idx_right);
        all_cue_idx = cat(2,all_cueleft_idx,all_cueright_idx);
        
end

%% debugging output
if debug_mode
    all_left_ts = EVTimeStamps(all_cueleft_idx);
    all_left_ts = all_left_ts*10^-6;
    
    all_right_ts = EVTimeStamps(all_cueright_idx);
    all_right_ts = all_right_ts*10^-6;
    
    figure;
    load(FindFile('*vt.mat'))
    plot(Data(x),Data(y),'.','MarkerSize',1,'Color',[0.5 0.5 0.5]); hold on;
    
    for iE = 1:length(all_left_ts)
        plot(Data(Restrict(x,all_left_ts(iE),all_left_ts(iE)+1)),Data(Restrict(y,all_left_ts(iE),all_left_ts(iE)+1)),'.r','MarkerSize',5,'Color',[1 0 0]);
    end
    
    for iE = 1:length(all_right_ts)
        plot(Data(Restrict(x,all_right_ts(iE),all_right_ts(iE)+1)),Data(Restrict(y,all_right_ts(iE),all_right_ts(iE)+1)),'.r','MarkerSize',5,'Color',[0 0 1]);
    end
    
    figure;
    plot(all_right_ts,0,'.r')
    hold on;
    plot(all_left_ts,0,'.b')
    
    % testing specific time intervals
    %     figure;
    %     plot(Data(x),Data(y),'.','MarkerSize',1,'Color',[0.5 0.5 0.5]); hold on;
    %     t_start = 1500; t_end = 1531;
    %
    %     for t0 = t_start:t_end-1
    %         t1 = t0+1;
    %         col = linspace(t0,t1,length(Data(Restrict(x,t0,t1)))); col = col-min(col); col = col./max(col);
    %         plot(Data(x),Data(y),'.','MarkerSize',1,'Color',[0.5 0.5 0.5]); hold on;
    %         try
    %             h = scatterplotC(Data(Restrict(x,t0,t1)),Data(Restrict(y,t0,t1)),col); %set(h,'MarkerSize',5);
    %         end
    %         pause; delete(h);
    %     end
    
end

%% collect everything
% we have:
%
%all_cueleft_idx 
%all_cueright_idx
%all_cue_idx
%
%c1_idx_left, etc...
%
% cue_correct_idx
% cue_correct_lat (for currect cues only)
%
% pb0_correct_idx and cue0_correct_idx (preceding cue)
% pb1_correct_idx and cue1_correct_idx (preceding cue)
%
% f0_idx
% f1_idx

% recall also that function input is:
%
%cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'
%cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
%cfg.trialdef.block = 'both'; % could be 'value', 'risk'
%cfg.trialdef.cue = {'c1','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'}

% first, get idx's of specified (eligible) cues
eligible_cue_idx = [];
for iCue = 1:length(cfg.trialdef.cue)
   
    cue_string = cat(2,cfg.trialdef.cue{iCue},'_idx');
    temp = eval(cue_string);
    eligible_cue_idx = cat(1,eligible_cue_idx,temp);
    
end

% exclude cues based on location
load(FindFile('*vt.mat'))
x_interp = interp1(Range(x),Data(x),EVTimeStamps(eligible_cue_idx)*10^-6);
keep_idx = find(x_interp > 150 & x_interp < 450);
if length(keep_idx) ~= length(eligible_cue_idx)
   fprintf('WARNING: %d cues removed based on location.\n',length(eligible_cue_idx)-length(keep_idx));
end
eligible_cue_idx = eligible_cue_idx(keep_idx);


switch cfg.trialdef.eventtype
    
    case 'nosepoke'
        
        switch cfg.trialdef.location
            case 'left'
                keep_idx = pb0_correct_idx; keep_cue_idx = cue0_correct_idx;
            case 'right'
                keep_idx = pb1_correct_idx; keep_cue_idx = cue1_correct_idx;
            case 'both'
                keep_idx = cat(2,pb0_correct_idx,pb1_correct_idx); 
                keep_cue_idx = cat(2,cue0_correct_idx,cue1_correct_idx); 
        end
        
        % now only include eligible cues
        [~,keep_a,~] = intersect(keep_cue_idx,eligible_cue_idx);
        keep_idx = keep_idx(keep_a); keep_cue_idx = keep_cue_idx(keep_a);
        
    case 'reward'
        
        switch cfg.trialdef.location
            case 'left'
                keep_idx = f0_idx; keep_cue_idx = cue0_correct_idx;
            case 'right'
                keep_idx = f1_idx; keep_cue_idx = cue1_correct_idx;
            case 'both'
                keep_idx = cat(1,f0_idx,f1_idx); 
                keep_cue_idx = cat(2,cue0_correct_idx,cue1_correct_idx); 
        end
        
        % now only include eligible cues
        [~,keep_a,~] = intersect(keep_cue_idx,eligible_cue_idx);
        keep_idx = keep_idx(keep_a); keep_cue_idx = keep_cue_idx(keep_a);
        
    case 'cue'
    
        switch cfg.trialdef.location
            case 'left'
                keep_idx = all_cueleft_idx;
            case 'right'
                keep_idx = all_cueright_idx;
            case 'both'
                keep_idx = all_cue_idx;
        end
        
        keep_idx = intersect(keep_idx,eligible_cue_idx);
        keep_idx = intersect(keep_idx,cue_correct_idx);
        
end

%% finally, generate output
keep_idx = keep_idx(~isnan(keep_idx)); % can be nan if e.g. a cue or feeder fire is somehow missing but lengths need to be maintained
event.idx = keep_idx;
switch OutputFormat
   
    case 'Timestamps'
        
        event.ts = EVTimeStamps(keep_idx);
        
    case 'Time'
    
        event.ts = EVTimeStamps(keep_idx) * 10^-6;
        
end


%%
if ~EvtOnly
    
    tvec = cfg.trialdef.hdr.FirstTimeStamp:cfg.trialdef.hdr.TimeStampPerSample:cfg.trialdef.hdr.LastTimeStamp;
    tvec = double(tvec);
    
    for iE = length(event.ts):-1:1
        
        event.idx(iE) = nearest(tvec,event.ts(iE));
        
    end
    
    % determine the number of samples before and after the trigger
    pretrig  = -round(cfg.trialdef.pre  * cfg.trialdef.hdr.Fs);
    posttrig =  round(cfg.trialdef.post * cfg.trialdef.hdr.Fs);
    
    % construct trl output
    trl = [];
    for j = 1:numel(event.idx)
        
        % find the sample corrsponding to this timestamp
        
        trlbegin = event.idx(j) + pretrig;
        trlend   = event.idx(j) + posttrig;
        offset   = pretrig;
        newtrl   = [trlbegin trlend offset];
        trl      = [trl; newtrl];
        
    end
    
else
    
    trl = [];
end


function [pb_out,f_out] = equalize_pbf(pb_in,f_in)
% add in missing feeder fires, can happen for certain sessions
next_pb = cat(2,pb_in(2:end),Inf);
for iPB = 1:length(pb_in)
    
   % check if a feeder fire comes before the next pb 
   next_f = find(f_in > pb_in(iPB) & f_in < next_pb(iPB));
   
   if ~isempty(next_f) % correct
       pb_out(iPB) = pb_in(iPB); f_out(iPB) = f_in(next_f(1));
   else % no feeder fire
       pb_out(iPB) = pb_in(iPB); f_out(iPB) = NaN;
   end
    
end
f_out = f_out';