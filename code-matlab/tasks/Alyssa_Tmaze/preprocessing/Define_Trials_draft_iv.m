% Define rest periods (pedestals)
% Define trial start-stop times


% written by AC. there's definitely a better way of doing the steps. and
% the variable naming is horrendous. should never be used for reference by
% other lab members (bad example, poor coding practice)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Get Session Epoch Intervals                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

% get ExpKeys info
LoadExpKeys

cfg.fc = ExpKeys.goodSWR(1);  

csc = LoadCSC(cfg);

tvec = csc.tvec;
data = csc.data;

pre_idx = nearest_idx3(ExpKeys.prerecord,tvec); % [start stop] indices for prerecord
task_idx = nearest_idx3(ExpKeys.task,tvec); % " for task
post_idx = nearest_idx3(ExpKeys.postrecord,tvec); % " for postrecord

% Grab all the indices for each epoch so they can be plotted
pre = pre_idx(1):pre_idx(2);
task = task_idx(1):task_idx(2);
post = post_idx(1):post_idx(2);

%% visually check the recording epochs

figure; hold on;

plot(tvec(pre),data(pre),'r',tvec(task),data(task),'g',tvec(post),data(post),'b');
legend('prerecord','main task','postrecord');

%% load vt
pos = LoadPos([]); 

% pull out x, y, timestamps:

xpos = getd(pos,'x'); % x coordinates
ypos = getd(pos,'y'); % y coordinates
tpos = pos.tvec; % timestamps

% cut position data based on epochs

% cut out prerecord and postrecord before working with main task
% pedestal-track interfaces

% but first, convert epoch times (taken from CSC) to pos times
pre_idx_pos = nearest_idx3(tvec(pre_idx),tpos)';
task_idx_pos = nearest_idx3(tvec(task_idx),tpos)';
post_idx_pos = nearest_idx3(tvec(post_idx),tpos)';

% now cut
% tidx = "time index" or "time indices"

pre_tidx = pre_idx_pos(1):pre_idx_pos(2);
task_tidx = task_idx_pos(1):task_idx_pos(2);
post_tidx = post_idx_pos(1):post_idx_pos(2);

%% Plot position data (visual check)

plot(xpos(pre_tidx),ypos(pre_tidx),'r');hold on; set(gca,'YDir','reverse');
plot(xpos(task_tidx),ypos(task_tidx),'g');
plot(xpos(post_tidx),ypos(post_tidx),'b');

%% plot 3D, incl time (visual check)
figure;
plot3(xpos(pre_tidx),ypos(pre_tidx),pre_tidx,'r');hold on;
plot3(xpos(task_tidx),ypos(task_tidx),task_tidx,'g');
plot3(xpos(post_tidx),ypos(post_tidx),post_tidx,'b');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Define rest periods (pedestals)                     %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% section is horribly written. inefficient and convoluted. don't use it for reference.  

%% Polygon alternative: left pedestal
xtask = xpos(task_idx_pos(1):task_idx_pos(2));
ytask = ypos(task_idx_pos(1):task_idx_pos(2));
plot(xtask,ytask,'k'); hold on; set(gca,'YDir','reverse');
maximize
title('Define n-point boundary containing left pedestal. Press Enter to commit.');
[xL,yL] = ginput;

% Polygon alternative: find if x and y points that are contained within the boundry (boolean)

onPedL = inpolygon(xtask,ytask,xL,yL); % logical. this checks to see if the points are contained within the polygon you drew
plot(xtask(onPedL),ytask(onPedL),'g');

% Polygon alternative: right pedestal
% draw a boundry for right ped
title('Define n-point boundary containing right pedestal. Press Enter to commit.');
[xR,yR] = ginput;
onPedR = inpolygon(xtask,ytask,xR,yR);
plot(xtask(onPedR),ytask(onPedR),'r');

%title('Define n-point boundary containing start box. Press Enter when done.');
%[xStart,yStart] = ginput;

%title('Define n-point boundary containing L stop box. Press Enter when done.');
%[xLstop,yLstop] = ginput;

%title('Define n-point boundary containing R stop box. Press Enter when done.');
%[xRstop,yRstop] = ginput;

title('Define n-point boundary containing track. Pay special attn to the start platform! Press Enter to commit.'); % track minus undesirable areas, like very start, plus rewards...on second thought, include start box b/c true beginning of trial
[xTrack,yTrack] = ginput;
% make sure the boundries do not overlap!
%inStartBox = inpolygon(xtask,ytask,xStart,yStart);
%inLstop = inpolygon(xtask,ytask,xLstop,yLstop);
%inRstop = inpolygon(xtask,ytask,xRstop,yRstop);
onPed = logical(onPedL + onPedR); % rat is on either pedestal. adding two logicals makes a double, so conv back so can be plotted.
onTrack = inpolygon(xtask,ytask,xTrack,yTrack);

%plot(xtask(inStartBox),ytask(inStartBox),'b');
%plot(xtask(inLstop),ytask(inLstop),'m');
%plot(xtask(inRstop),ytask(inRstop),'c');
plot(xtask(onTrack),ytask(onTrack),'b');
%plot(xtask(onPed),ytask(onPed),'y');
title('Run cell again to redo');

clear xStart yStart xLstop yLstop xRstop yRstop xL yL xR yR


%% To view above final plot again without having to re-run the cell:
figure; hold on;
plot(xtask,ytask,'k'); set(gca,'YDir','reverse');
plot(xtask(onPedL),ytask(onPedL),'g');
plot(xtask(onPedR),ytask(onPedR),'r');
plot(xtask(onTrack),ytask(onTrack),'b');
title('Task boundaries');

%% get T-maze photobeam timestamps and reward delivery timestamps
% due to messy boundaries, there might be many on-off times that aren't
% real. For this reason, we also need to use the reward firing times so
% we have some way of finding the first entry or exit during a single rest period

pb_evt = getEvents_Tmaze(); % photobeam events (and reward delivery)
ts_L = pb_evt.food_dispensed; % timestamps for L reward
ts_R = pb_evt.water_dispensed; % timestamps for R reward
%ts_pbL = pb_evt.left_pb; % timestamps for L photobeam
ts_pbC = pb_evt.center_pb; % '' center photobeam
%ts_pbR = pb_evt.right_pb; % '' R photobeam


% condense ts down to single passage events, or w/e
% find where the difference in timestamps is larger than 1s. We'll say this
% border contitues different passages. However, this misses the last
% passage, so keep the last timestamp to represent the last passage. Prob
% won't make any sense to me a few months from now...
%idx_pbL = cat(2,find(diff(pbLts)>1),length(pbLts));
idx_pbC = cat(2,find(diff(ts_pbC)>1),length(ts_pbC));
%idx_pbR = cat(2,find(diff(pbRts)>1),length(pbRts));
% the above arrays contain indices of passage times

% don't know what to call this next bit:
% make a matrix containing zeros in place of ts. replace the zeros with 1s
% where L pb is broken, 2s for C pb, and 3s for R pb.

tvec_task = tvec(task); % restrict tvec to task times only

checkpoint_id = zeros(size(tvec_task));


checkpoint_id(nearest_idx3(ts_L,tvec_task)) = 1; % convert values to 1 where pbL is activated
checkpoint_id(nearest_idx3(ts_pbC(idx_pbC),tvec_task)) = 2; % convert values to 2 where pbC is activated
checkpoint_id(nearest_idx3(ts_R,tvec_task)) = 3; % convert values to 3 where pbR is activated
% the above generates an array with the identity of each checkpoint preserved
% in order. It's convoluted but the only way I could understand to do this. 

%% check checkpoint_id (should see alternating 2,1,2,3 and so on, corresponding to choices)
figure; hold on; set(gca,'YLim',[-0.5 3.5]); 
plot(checkpoint_id); 
 
%% try with checkpoint_id? get rest period boundaries (on pedestals) in terms of ts
% can also use this to generate trial sequence
checkpoint_idx = find(checkpoint_id); % all the indices where a rat passed a certain checkpoint
bound_Lped = []; % to hold the Lped time interval boundaries
bound_Rped = []; % to hold the Rped time interval boundaries
bound_allPed = [];
sequence = {};

for iCheckpoint = 1:length(checkpoint_idx)-1
    
    
        checkpointA = checkpoint_id(checkpoint_idx(iCheckpoint));
        checkpointB = checkpoint_id(checkpoint_idx(iCheckpoint+1));
    
        if checkpointA == 1 && checkpointB == 2 % it's a left-ped boundary (because 1 reprs L reward and 2 reprs pbC)
        
            bound_Lped = cat(1,bound_Lped,[tvec_task(checkpoint_idx(iCheckpoint)) tvec_task(checkpoint_idx(iCheckpoint+1))]);
            sequence = cat(2,sequence,'L');
        end
    
        if checkpointA == 3 && checkpointB == 2 %% it's a right-ped boundary (because 3 reprs R reward and 2 reprs pbC)
        
            bound_Rped = cat(1,bound_Rped,[tvec_task(checkpoint_idx(iCheckpoint)) tvec_task(checkpoint_idx(iCheckpoint+1))]);
            sequence = cat(2,sequence,'R');
        end
        
        % handle endpoint case
        if iCheckpoint == length(checkpoint_idx) - 1
            checkpointA = checkpoint_id(checkpoint_idx(iCheckpoint + 1)); 
            
            if checkpointA == 1 % Left ped finish
                bound_Lped = cat(1,bound_Lped,[tvec_task(checkpoint_idx(iCheckpoint+1)) tvec_task(length(tvec_task))]);
                sequence = cat(2,sequence,'L');
            end
            
            if checkpointA == 3 % right ped finish
            bound_Rped = cat(1,bound_Rped,[tvec_task(checkpoint_idx(iCheckpoint+1)) tvec_task(length(tvec_task))]);
            sequence = cat(2,sequence,'R');
            end
        end       
   
end

% convert to VT time intervals
% reminder: these intervals are between feeder fires and center pb breaks.
% they are not the actual times the rat was on the ped, but these are
% contained within this interval. 

ttask = tpos(task_tidx); % timestamps for task

tvt_onPedL = [nearest_idx3(bound_Lped(:,1),ttask) nearest_idx3(bound_Lped(:,2),ttask)];
tvt_onPedR = [nearest_idx3(bound_Rped(:,1),ttask) nearest_idx3(bound_Rped(:,2),ttask)];
%tvt_onPed = [nearest_idx3(bound_allPed(:,1),ttask) nearest_idx3(bound_allPed(:,2),ttask)];

%tvt_onPedL = [tpos(nearest_idx3(bound_Lped(:,1),tpos))' tpos(nearest_idx3(bound_Lped(:,2),tpos))']; %timestamps
%tvt_onPedR = [tpos(nearest_idx3(bound_Rped(:,1),tpos))' tpos(nearest_idx3(bound_Rped(:,2),tpos))']; % timestamps

% get the actual onPed intervals

% for pedL

ttask_onPedL = ttask(onPedL);
tidx_onPedL = [];

for itvt = 1:length(tvt_onPedL) % itvt will be the row we're working with
    iv_temp = tvt_onPedL(itvt,:); % take all of the itvt'th row
    keep = [nearval3(ttask(iv_temp(1)),ttask_onPedL,1) nearval3(ttask(iv_temp(2)),ttask_onPedL,-1)]; % interval for keeping pos data 
    % direction in nearest is important, because the rat takes longer to get on the ped than he does to get from center photobeam to track arm
    keep = [find(ttask == keep(1)) find(ttask == keep(2))]; % conv back to idx
    tidx_onPedL = cat(1,tidx_onPedL,keep); % this in terms of position data, rather than nlx data         
end

% for pedR
ttask_onPedR = ttask(onPedR);
tidx_onPedR = [];

for itvt = 1:length(tvt_onPedR) % itvt will be the row we're working with
    iv_temp = tvt_onPedR(itvt,:); % take all of the itvt'th row
    keep = [nearval3(ttask(iv_temp(1)),ttask_onPedR,1) nearval3(ttask(iv_temp(2)),ttask_onPedR,-1)]; % interval for keeping pos data
    keep = [find(ttask == keep(1)) find(ttask == keep(2))]; % conv back to idx
    tidx_onPedR = cat(1,tidx_onPedR,keep); % this in terms of position data, rather than nlx data         
end

% omfg what a pain this whole section was


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   Get trial intervals                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ttask_track = ttask(onTrack); % timestamps on track

% make a single matrix with trial start and stop times based on a
% cat of the L and R ped times. pair trial iv with sequence to pull out
% which is L vs R.

tidx_onPed = sortrows(cat(1,tidx_onPedL,tidx_onPedR)); % sortrows: mon sauveur?

% handle first trial 
keep = [nearval3(ttask(1),ttask_track) nearval3(ttask(tidx_onPed(1)),ttask_track,-1)];
%bound(1) = find(ttask == bound(1));
%bound(2) = nearval3(ttask(tidx_onPed(1)),ttask_track,1);
%bound(2) = find(ttask == bound);
keep = [find(ttask == keep(1)) find(ttask == keep(2))];
tidx_onTrack = [keep(1) keep(2)];

% handle additional trials
for irow = 1:length(tidx_onPed)-1 % #trials = # rows
    bound = [ttask(tidx_onPed(irow)) ttask(tidx_onPed(irow+1))]; % These will be the boundaries for pulling out pos data per trial
    keep = [nearval3(bound(1),ttask_track,1) nearval3(bound(2),ttask_track,-1)];
    keep = [find(ttask == keep(1)) find(ttask == keep(2))]; % conv back to idx
    tidx_onTrack = cat(1,tidx_onTrack,keep);
end

% separate into L and R trials
trial_iv_idx_L = [];
trial_iv_idx_R = [];
for iTrial = 1:length(sequence)
    if strcmp(sequence(iTrial),'L')
        trial_iv_idx_L = cat(1,trial_iv_idx_L,tidx_onTrack(iTrial,:));
    end
    if strcmp(sequence(iTrial),'R')
        trial_iv_idx_R = cat(1,trial_iv_idx_R,tidx_onTrack(iTrial,:));
    end
end

%% PLOT THE TRIALS

%tiv_Ltrials = []; % to hold the Ltrials time interval
%tiv_Rtrials = []; % to hold the Rtrials time interval

% invent a convoluted way to figure out which trials are forced or which
% need to be discarded
forcedTrials = ExpKeys.forcedTrials;
discardTrials = ExpKeys.badTrials; 
trialType = zeros(size(sequence)); % ID array. contains 1 for forced, and 2 for discard
trialType(forcedTrials) = 1;
trialType(discardTrials) = 2;

% set number of plots
numRow = 4; numCol = 5;
if length(tidx_onTrack)>20 % if we did more than 20 trials, use these parameters instead
    numRow = 5;
    numCol = 5;
end

% plot the trials
figure; %hold on;

for iTrial = 1:length(tidx_onTrack)
    iv_temp = tidx_onTrack(iTrial,:);
    subplot(numRow,numCol,iTrial); hold on;
    plot(xtask,ytask,'k');set(gca,'YDir','reverse');
    color = 'g'; % plot on green for normal trials
    type = ''; % empty string for normal trials
    if trialType(iTrial) == 1 % set variables for plotting forced trials
        color = 'c';
        type = 'forced';
    end
    
    if trialType(iTrial) == 2 %set variables for plotting discard trials
        color = 'r';
        type = 'discard';
    end
    plot(xtask(iv_temp(1):iv_temp(2)),ytask(iv_temp(1):iv_temp(2)),color,'LineWidth',3);
    dir = sequence{iTrial}; % direction (L or R); need to use {} instead of () or plot title is weird
    title_str = ['Trial ',num2str(iTrial),': ',dir, ' ',type];
    title(title_str);
end
maximize

%% plot the rest periods
figure; 

% set number of plots
numRow = 4; numCol = 5;
if length(sequence)>20 % if we did more than 20 trials, use these parameters instead
    numRow = 5;
    numCol = 5;
end

for iRest = 1:length(tidx_onPed)
    iv_temp = tidx_onPed(iRest,:);
    subplot(numRow,numCol,iRest); hold on;
    plot(xtask,ytask,'k');set(gca,'YDir','reverse');
    dur = (ttask(iv_temp(2))-ttask(iv_temp(1)))/60;
    plot(xtask(iv_temp(1):iv_temp(2)),ytask(iv_temp(1):iv_temp(2)),'m','LineWidth',3);
    dir = sequence{iRest}; % direction (L or R); need to use {} instead of () or plot title is weird
    title_str = ['Rest ',num2str(iRest),': ',dir,', ',num2str(dur),' min'];
    title(title_str);
end
maximize

%% if you're happy with the results, i.e.:

% THE TRIALS SHOULD NOT INCLUDE ANY PART OF THE PEDESTALS; IF THEY DO YOU
% DREW THE TRACK BOUNDARY INCORRECLTY. MAKE SURE YOU DRAW THE OUTLINE OF THE
% STARTING PLATFORM, WHICH IS SLIGHTLY WIDER THAN THE MAIN TRACK. SEE LAB
% NOTES FOR A REMINDER.

% trial intervals:
trial_iv_ts = ttask(tidx_onTrack);
trial_iv_ts = iv(trial_iv_ts(:,1),trial_iv_ts(:,2));

% left trials
trial_iv_L = ttask(trial_iv_idx_L);
trial_iv_L = iv(trial_iv_L(:,1),trial_iv_L(:,2));

% right trials
trial_iv_R = ttask(trial_iv_idx_R);
trial_iv_R = iv(trial_iv_R(:,1),trial_iv_R(:,2));

% rest intervals for both pedestals together:
rest_iv_ts = ttask(tidx_onPed);
rest_iv_ts = iv(rest_iv_ts(:,1),rest_iv_ts(:,2));

% rest intervals for left pedestal:
rest_iv_ts_pedL = ttask(tidx_onPedL);
rest_iv_ts_pedL = iv(rest_iv_ts_pedL(:,1),rest_iv_ts_pedL(:,2));

% rest intervals for right pedestal:
rest_iv_ts_pedR = ttask(tidx_onPedR);
rest_iv_ts_pedR = iv(rest_iv_ts_pedR(:,1),rest_iv_ts_pedR(:,2));

%% save as a metadata field called "taskvars"

% WARNING: running this section overwrites existing taskvars field (if one exists already)!

taskvars = struct('trial_iv',trial_iv_ts,'trial_iv_L',trial_iv_L,'trial_iv_R',trial_iv_R,'rest_iv',rest_iv_ts,'rest_iv_pedL',rest_iv_ts_pedL,'rest_iv_pedR',rest_iv_ts_pedR);
taskvars.sequence = sequence;

% first check if metadata exists yet
loaded = LoadMetadata2;

if ~loaded % then it doesn't exist yet, so make a metadata struct
    metadata.taskvars = taskvars;
else % it does, so add a new field
    metadata.taskvars = taskvars;
end

[~,name,~] = fileparts(pwd); % pwd is your current folder, we just want its namepart

% now save
savename = strcat(name,'-metadata.mat'); % use the folder's name, but concatenate it with '-metadata'
save(savename,'metadata'); % this saves the specified variables under the given [save]name