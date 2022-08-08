% Script to generate sinusoidal waves
% TODO: Generalize to do this for various step sizes and input frequencies
clear;
close all;
n_phases = 4;
c_ord = linspecer(n_phases);
f1 = 100;
cycles = 8; % Total cycles of sinewave
x = -cycles*pi:pi/(f1*n_phases):cycles*pi;
step_size = floor(length(x)/8);
%%
wave_phase_1 = sin(x);
slope_1 = diff(wave_phase_1);
wave_phase_2 = sin(1.5*x);
slope_2 = diff(wave_phase_2);
wave_phase_3 = sin(1.8*x);
slope_3 = diff(wave_phase_3);
first_wave = zeros(size(x));
first_wave(1:step_size) = wave_phase_1(1:step_size);
search_among = wave_phase_2;
if slope_1(step_size) >= 0
    search_among(slope_2<0) = nan;
else
    search_among(slope_2>=0) = nan;
end
next_idx =  nearest_idx2(first_wave(step_size), search_among(1:end-step_size));

first_wave(step_size+1:2*step_size) = wave_phase_2(next_idx+1:next_idx+step_size);
search_among = wave_phase_3;
if slope_2(next_idx+step_size) >= 0
    search_among(slope_3<0) = nan;
else
    search_among(slope_3>=0) = nan;
end
next_idx =  nearest_idx2(first_wave(2*step_size), search_among(1:end-step_size));

first_wave(2*step_size+1:3*step_size) = wave_phase_3(next_idx+1:next_idx+step_size);
search_among = wave_phase_1;
if slope_3(next_idx+step_size) >= 0
    search_among(slope_1<0) = nan;
else
    search_among(slope_1>=0) = nan;
end
next_idx =  nearest_idx2(first_wave(3*step_size), search_among(1:end-step_size));

first_wave(3*step_size+1:4*step_size) = wave_phase_1(next_idx+1:next_idx+step_size);
search_among = wave_phase_2;
if slope_1(next_idx+step_size) >= 0
    search_among(slope_2<0) = nan;
else
    search_among(slope_2>=0) = nan;
end
next_idx =  nearest_idx2(first_wave(4*step_size), search_among(1:end-step_size));

first_wave(4*step_size+1:5*step_size) = wave_phase_2(next_idx+1:next_idx+step_size);
search_among = wave_phase_3;
if slope_2(next_idx+step_size) >= 0
    search_among(slope_3<0) = nan;
else
    search_among(slope_3>=0) = nan;
end
next_idx =  nearest_idx2(first_wave(5*step_size), search_among(1:end-step_size));

first_wave(5*step_size+1:6*step_size) = wave_phase_3(next_idx+1:next_idx+step_size);
search_among = wave_phase_1;
if slope_3(next_idx+step_size) >= 0
    search_among(slope_1<0) = nan;
else
    search_among(slope_1>=0) = nan;
end
next_idx =  nearest_idx2(first_wave(6*step_size), search_among(1:end-step_size));

first_wave(6*step_size+1:7*step_size) = wave_phase_1(next_idx+1:next_idx+step_size);
search_among = wave_phase_2;
if slope_1(next_idx+step_size) >= 0
    search_among(slope_2<0) = nan;
else
    search_among(slope_2>=0) = nan;
end
next_idx =  nearest_idx2(first_wave(7*step_size), search_among(1:end-step_size));

first_wave(7*step_size+1:end) = wave_phase_2(next_idx+1:next_idx+step_size+1);

%%
color_scheme = zeros(size(x));
phase_angles = angle(hilbert(first_wave));
for i = 1:length(color_scheme)
    if phase_angles(i) <= -0.5*pi
        color_scheme(i) = 1;
    elseif phase_angles(i) <= 0
        color_scheme(i) = 2;
    elseif phase_angles(i) <= 0.5*pi
        color_scheme(i) = 3;
    else
        color_scheme(i) = 4;           
    end
end
changes_in_phase = find(diff(color_scheme));
figure;
subplot(2,1,1)
hold on;
j = 1;
for i=1:length(changes_in_phase)
    plot(j:changes_in_phase(i), first_wave(j:changes_in_phase(i)),...
        'color', c_ord(color_scheme(changes_in_phase(i)),:), 'linewidth', 5);
    j = changes_in_phase(i)+1;
end
if j<length(color_scheme)
    plot(j:length(color_scheme), first_wave(j:end), ...
        'color', c_ord(color_scheme(end),:), 'linewidth', 5);
end
%%
second_wave = sin(x);
color_scheme = zeros(size(x));
phase_angles = angle(hilbert(second_wave));
for i = 1:length(color_scheme)
    if phase_angles(i) <= -0.5*pi
        color_scheme(i) = 1;
    elseif phase_angles(i) <= 0
        color_scheme(i) = 2;
    elseif phase_angles(i) <= 0.5*pi
        color_scheme(i) = 3;
    else
        color_scheme(i) = 4;           
    end
end
changes_in_phase = find(diff(color_scheme));
j = 1;
for i=1:length(changes_in_phase)
    plot(j:changes_in_phase(i), second_wave(j:changes_in_phase(i))-2.2,...
        'color', c_ord(color_scheme(changes_in_phase(i)),:), 'linewidth', 5);
    j = changes_in_phase(i)+1;
end
if j<length(color_scheme)
    plot(j:length(color_scheme), second_wave(j:end)-2.2, ...
        'color', c_ord(color_scheme(end),:), 'linewidth', 5);
end


%%
wave_phase_1 = sin((pi/4)+x);
slope_1 = diff(wave_phase_1);
wave_phase_2 = sin(1.5*x);
slope_2 = diff(wave_phase_2);
wave_phase_3 = sin((pi/4)+x);
slope_3 = diff(wave_phase_3);

step_size = 2*step_size;
third_wave = zeros(size(x));
third_wave(1:step_size) = wave_phase_1(1:step_size);
search_among = wave_phase_2;
if slope_1(step_size) >= 0
    search_among(slope_2<0) = nan;
else
    search_among(slope_2>=0) = nan;
end
next_idx =  nearest_idx2(third_wave(step_size), search_among(1:end-step_size));

third_wave(step_size+1:2*step_size) = wave_phase_2(next_idx+1:next_idx+step_size);
search_among = wave_phase_3;
if slope_2(next_idx+step_size) >= 0
    search_among(slope_3<0) = nan;
else
    search_among(slope_3>=0) = nan;
end
next_idx =  nearest_idx2(third_wave(2*step_size), search_among(1:end-step_size));

third_wave(2*step_size+1:3*step_size) = wave_phase_3(next_idx+1:next_idx+step_size);
search_among = wave_phase_1;
if slope_3(next_idx+step_size) >= 0
    search_among(slope_1<0) = nan;
else
    search_among(slope_1>=0) = nan;
end
next_idx =  nearest_idx2(third_wave(3*step_size), search_among(1:end-step_size));

third_wave(3*step_size+1:end) = wave_phase_1(next_idx+1:next_idx+step_size+1);

%%
color_scheme = zeros(size(x));
phase_angles = angle(hilbert(third_wave));
for i = 1:length(color_scheme)
    if phase_angles(i) <= -0.5*pi
        color_scheme(i) = 1;
    elseif phase_angles(i) <= 0
        color_scheme(i) = 2;
    elseif phase_angles(i) <= 0.5*pi
        color_scheme(i) = 3;
    else
        color_scheme(i) = 4;           
    end
end
changes_in_phase = find(diff(color_scheme));
j = 1;
for i=1:length(changes_in_phase)
    plot(j:changes_in_phase(i), third_wave(j:changes_in_phase(i))-4.5,...
        'color', c_ord(color_scheme(changes_in_phase(i)),:), 'linewidth', 5);
    j = changes_in_phase(i)+1;
end
if j<length(color_scheme)
    plot(j:length(color_scheme), third_wave(j:end)-4.5, ...
        'color', c_ord(color_scheme(end),:), 'linewidth', 5);
end
