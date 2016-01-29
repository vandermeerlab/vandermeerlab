function [ times ] = phase_times(csc_file, rat_id)
% Takes in a csc file as input and determines where the breaks in the
% experiment occured for different phases. Returns a times object that
% contains each phase with a start(1) and end(2) time (in seconds). 
% - Requires emi_loadCSC.m
% prerecord: rat on pedestal, usually around 5min (300sec)
% phase1: "U"-shaped track, usually around 8min (480sec)
% pauseA: rat on pedestal, set-up shortcut, usually around 10min (600sec)
% phase2: barriers, usually around 20min (1200sec)
% pauseB: rat on pedestal, barriers removed at the beginning, usually around 30min (1800sec)
% phase3: barriers, usually around 50min (3000sec)
% postrecord: rat on pedestal, usually around 5min (300sec)
% - Returns times of all these experimental phases.
%
% Usage:
% rat_id = 'R068_EI';
% unique_folder = 'R068-2014-12-01_recording';
% cd(fullfile('C:\Users\Emily\Desktop\',rat_id,unique_folder));
% times = phase_times('R068-2014-12-01-CSC10d.ncs', 'R066_day1');

% Load CSC from csc_file
cfg = []; cfg.fc = {csc_file};
csc = LoadCSC(cfg);

% Determine when gap between two time points is 10x longer than usual
dt = csc.tvec(2) - csc.tvec(1);
time = csc.tvec(2:end) - csc.tvec(1:end-1);
pause = find(time > dt*10);


if length(pause) == 6
    idx_start_prerecord = 1;
    idx_end_prerecord = pause(1);
    idx_start_phase1 = pause(1)+1;
    idx_end_phase1 = pause(2);
    idx_start_pause_a = pause(2)+1;
    idx_end_pause_a = pause(3);
    idx_start_phase2 = pause(3)+1;
    idx_end_phase2 = pause(4);
    idx_start_pause_b = pause(4)+1;
    idx_end_pause_b = pause(5);
    idx_start_phase3 = pause(5)+1;
    idx_end_phase3 = pause(6);
    idx_start_postrecord = pause(6)+1;
    idx_end_postrecord = max(size(csc.tvec));
elseif length(pause) == 8
    idx_start_one = 1;
    idx_end_one = pause(1);
    idx_start_two = pause(1)+1;
    idx_end_two = pause(2);
    idx_start_three = pause(2)+1;
    idx_end_three = pause(3);
    idx_start_four = pause(3)+1;
    idx_end_four = pause(4);
    idx_start_five = pause(4)+1;
    idx_end_five = pause(5);
    idx_start_six = pause(5)+1;
    idx_end_six = pause(6);
    idx_start_seven = pause(6)+1;
    idx_end_seven = max(7);
    idx_start_eight = pause(7)+1;
    idx_end_eight = pause(8);
    idx_start_nine = pause(8)+1;
    idx_end_nine = max(size(csc.tvec));
elseif length(pause) == 5
    idx_start_one = 1;
    idx_end_one = pause(1);
    idx_start_two = pause(1)+1;
    idx_end_two = pause(2);
    idx_start_three = pause(2)+1;
    idx_end_three = pause(3);
    idx_start_four = pause(3)+1;
    idx_end_four = pause(4);
    idx_start_five = pause(4)+1;
    idx_end_five = pause(5);
    idx_start_six = pause(5)+1;
    idx_end_six = max(size(csc.tvec));
elseif length(pause) == 9
    idx_start_one = 1;
    idx_end_one = pause(1);
    idx_start_two = pause(1)+1;
    idx_end_two = pause(2);
    idx_start_three = pause(2)+1;
    idx_end_three = pause(3);
    idx_start_four = pause(3)+1;
    idx_end_four = pause(4);
    idx_start_five = pause(4)+1;
    idx_end_five = pause(5);
    idx_start_six = pause(5)+1;
    idx_end_six = pause(6);
    idx_start_seven = pause(6)+1;
    idx_end_seven = pause(7);
    idx_start_eight = pause(7)+1;
    idx_end_eight = pause(8);
    idx_start_nine = pause(8)+1;
    idx_end_nine = pause(9);
    idx_start_ten = pause(9)+1;
    idx_end_ten = max(size(csc.tvec));
else
    print('check recording times')
end


% Store start and end times in times object
if length(pause) == 6
    times.prerecord = [csc.tvec(idx_start_prerecord), csc.tvec(idx_end_prerecord)];
    times.phase1 = [csc.tvec(idx_start_phase1), csc.tvec(idx_end_phase1)];
    times.pauseA = [csc.tvec(idx_start_pause_a), csc.tvec(idx_end_pause_a)];
    times.phase2 = [csc.tvec(idx_start_phase2), csc.tvec(idx_end_phase2)];
    times.pauseB = [csc.tvec(idx_start_pause_b), csc.tvec(idx_end_pause_b)];
    times.phase3 = [csc.tvec(idx_start_phase3), csc.tvec(idx_end_phase3)];
    times.postrecord = [csc.tvec(idx_start_postrecord), csc.tvec(idx_end_postrecord)];

    % Store idx
    times.idx.prerecord = [idx_start_prerecord,idx_end_prerecord];
    times.idx.phase1 = [idx_start_phase1,idx_end_phase1];
    times.idx.pauseA = [idx_start_pause_a,idx_end_pause_a];
    times.idx.phase2 = [idx_start_phase2,idx_end_phase2];
    times.idx.pauseB = [idx_start_pause_b,idx_end_pause_b];
    times.idx.phase3 = [idx_start_phase3,idx_end_phase3];
    times.idx.postrecord = [idx_start_postrecord,idx_end_postrecord];
end

if rat_id == 'R067_day1'
    times.prerecord = [csc.tvec(idx_start_one), csc.tvec(idx_end_one)];
    times.phase1 = [csc.tvec(idx_start_two), csc.tvec(idx_end_two)];
    times.pauseA = [csc.tvec(idx_start_three), csc.tvec(idx_end_three)];
    times.phase2 = [csc.tvec(idx_start_four), csc.tvec(idx_end_four)];
    times.pauseB = [csc.tvec(idx_start_five), csc.tvec(idx_end_five)];
    times.phase3 = [csc.tvec(idx_start_six), csc.tvec(idx_end_eight)];
    times.postrecord = [csc.tvec(idx_start_nine), csc.tvec(idx_end_nine)];
end

if rat_id == 'R067_day6'
    times.prerecord = [csc.tvec(idx_start_one), csc.tvec(idx_end_one)];
    times.phase1 = [csc.tvec(idx_start_two), csc.tvec(idx_end_two)];
    times.pauseA = [csc.tvec(idx_start_three), csc.tvec(idx_end_three)];
    times.phase2 = [csc.tvec(idx_start_four), csc.tvec(idx_end_four)-1.87999e+03];
    times.pauseB = [csc.tvec(idx_end_four)-1.88e+03, csc.tvec(idx_end_four)];
    times.phase3 = [csc.tvec(idx_start_five), csc.tvec(idx_end_five)];
    times.postrecord = [csc.tvec(idx_start_six), csc.tvec(idx_end_six)];
end


if rat_id == 'R066_day1'
    times.prerecord = [csc.tvec(idx_start_one), csc.tvec(idx_end_one)];
    times.phase1 = [csc.tvec(idx_start_two), csc.tvec(idx_end_two)];
    times.pauseA = [csc.tvec(idx_start_three), csc.tvec(idx_end_four)];
    times.phase2 = [csc.tvec(idx_start_seven), csc.tvec(idx_end_seven)];
    times.pauseB = [csc.tvec(idx_start_eight), csc.tvec(idx_end_eight)];
    times.phase3 = [csc.tvec(idx_start_nine), csc.tvec(idx_end_nine)];
    times.postrecord = [csc.tvec(idx_start_ten), csc.tvec(idx_end_ten)];
end
  