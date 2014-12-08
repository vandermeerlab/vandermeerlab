%%%
%
% This script removes the reward/chewing/licking epoch from tetrodes
% needed for R042's data because was recorded relative to ground
%
% MvdM 13
%
%%%

%% parameters
pretrial_offset = -2; % start trial relative to last center pb break before reward
posttrial_offset = 0.5;
prepedestal_offset = 2;
postpedestal_offset = -2;

%% first get events
evt = getEvents_Tmaze();

%% load vt
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT('Vt1.nvt', [1 1 1 1 1 1], 1, 1, [] );
X(X == 0) = NaN; Y(Y == 0) = NaN; % this is useful so that (0,0) doesn't mess up position-based trial events later

%% user selected input for pre and postrecord
plot(Timestamps,X,'.','Color',[0.5 0.5 0.5],'MarkerSize',1);
hold on;

ref_axes = axis;

title('zoom in for prerecord start -- press a key to continue');
pause;
prerecord_start = ginput(1); prerecord_start = prerecord_start(1);

pre_h = plot([prerecord_start prerecord_start],[min(X) max(X)],'k');
axis(ref_axes);

title('zoom in for prerecord end -- press a key to continue');
pause;
prerecord_end = ginput(1); prerecord_end = prerecord_end(1);

pre_h = plot([prerecord_end prerecord_end],[min(X) max(X)],'k');
axis(ref_axes);

title('zoom in for postrecord start -- press a key to continue');
pause;
postrecord_start = ginput(1); postrecord_start = postrecord_start(1);

pre_h = plot([postrecord_start postrecord_start],[min(X) max(X)],'k');
axis(ref_axes);

title('zoom in for postrecord end -- press a key to continue');
pause;
postrecord_end = ginput(1); postrecord_end = postrecord_end(1);

pre_h = plot([postrecord_end postrecord_end],[min(X) max(X)],'k');
axis(ref_axes); title('');

title('press key to continue');
pause;

%% user selected input for pedestal
clf;

plot(X,Y,'.','Color',[0.5 0.5 0.5],'MarkerSize',1);
hold on;

title('define X coordinate for pedestal entered');
pedestal_enterX = ginput(1); pedestal_enterX  = pedestal_enterX(1);
pre_h = plot([pedestal_enterX pedestal_enterX],[min(Y) max(Y)],'k');

title('define X coordinate for pedestal exit');
pedestal_exitX = ginput(1); pedestal_exitX  = pedestal_exitX(1);
pre_h = plot([pedestal_exitX pedestal_exitX],[min(Y) max(Y)],'k');

title('press key to continue');
pause; close all;

%% get trial start and end times
run_end = cat(2,evt.water_dispensed,evt.food_dispensed); % note these are in seconds
run_end = sort(run_end,'ascend')+posttrial_offset;

figure;
clear run_start;
for iRun = 1:length(run_end)
   
    preceding_starts = evt.center_pb(evt.center_pb < run_end(iRun));
    run_start(iRun) = max(preceding_starts) + pretrial_offset;
    
    subplot(5,5,iRun);
    
    plot(X,Y,'.','Color',[0.5 0.5 0.5],'MarkerSize',1);
    hold on;
    
    trial_idx = find(Timestamps > run_start(iRun)*10^6 & Timestamps < run_end(iRun)*10^6);
    plot(X(trial_idx),Y(trial_idx),'.r','MarkerSize',10);
    title(sprintf('trial %d',iRun));
    
    axis off;
    
end

title('press key to continue');
pause; close all;

%% convert pedestal coordinates to enter and exit times

% first find all pedestal enters and exits
pedestal_entered = (X < pedestal_enterX);
pedestal_entered = diff(pedestal_entered); % transition from 0 to 1 means go from right to left onto pedestal
pedestal_entered = Timestamps(pedestal_entered == 1); % may be one sample off but we don't care

pedestal_exited = (X < pedestal_exitX);
pedestal_exited = diff(pedestal_exited); % transition from 0 to 1 means go from right to left onto pedestal
pedestal_exited = Timestamps(pedestal_exited == 1); % may be one sample off but we don't care

figure;
clear ped_start ped_end;
for iRun = 1:length(run_end)-1
   
    following_enters = pedestal_entered(pedestal_entered > run_end(iRun)*10^6);
    ped_start(iRun) = min(following_enters) + prepedestal_offset*10^6;
    
    following_exits = pedestal_exited(pedestal_exited > run_end(iRun)*10^6);
    ped_end(iRun) = min(following_exits) + postpedestal_offset*10^6;
    
    dt = ped_end(iRun)*10^-6 - ped_start(iRun)*10^-6;
    
    subplot(5,5,iRun);
    
    plot(X,Y,'.','Color',[0.5 0.5 0.5],'MarkerSize',1);
    hold on;
    
    trial_idx = find(Timestamps > ped_start(iRun) & Timestamps < ped_end(iRun));
    plot(X(trial_idx),Y(trial_idx),'.r','MarkerSize',10);
    title(sprintf('t %d (%.1f s)',iRun,dt));
    
    axis off;
    
end

issues = 0;
if issues
% ped_exclude = [6 8 11];
   for iE = 1:length(ped_exclude), ped_end(ped_exclude(iE)) = ped_start(ped_exclude(iE)); end
   
end

%% compile times
t_start = []; t_end = [];

t_start = cat(2,t_start,prerecord_start,postrecord_start); % in Timestamps
t_end = cat(2,t_end,prerecord_end,postrecord_end);

t_start = cat(2,t_start,run_start*10^6); % in seconds!
t_end = cat(2,t_end,run_end*10^6);

t_start = cat(2,t_start,ped_start);
t_end = cat(2,t_end,ped_end);

if any((t_start < t_end) == 0)
   disp('WARNING, times don''t match!'); 

else
   disp('Times seem ok');
   
   t_start = sort(t_start,'ascend');
   t_end = sort(t_end,'ascend');
   
   save('times.mat','t_start','t_end','prerecord_start','prerecord_end','postrecord_start','postrecord_end','run_start','run_end','ped_start','ped_end');
end


%% load test tt
fd = FindFiles('*.ntt');
fixTT12 = 0;

for iF = 1:length(fd)
    
    fname = fd{iF};
    
    disp(sprintf('Processing tt %s (%d/%d)...',fname,iF,length(fd)));
    
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(fname, [1 1 1 1 1], 1, 1, [] );
    
    if iF == 11 & fixTT12
        disp('***FIXING...');
        Samples(:,1,:) = Samples(:,1,:).*2; % if gain was set incorrectly
    end
    
    %% restrict
    keep_idx = [];
    for iR = 1:length(t_start)
        
        keep_idx = cat(2,keep_idx,find(Timestamps > t_start(iR) & Timestamps < t_end(iR)));
        
    end
    
    %% export
    fname_out = regexprep(fname,'\.','r\.');
    Mat2NlxSpike(fname_out, 0, 1, [], [1 1 1 1 1], Timestamps(keep_idx), ScNumbers(keep_idx), CellNumbers(keep_idx), Features(:,keep_idx), Samples(:,:,keep_idx), Header);
    
end