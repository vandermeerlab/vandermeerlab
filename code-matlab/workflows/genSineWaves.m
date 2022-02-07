close all;
n_phases = 5;
c_ord = linspecer(n_phases);
figure
f1 = 100;
cycles = 8;
x = -cycles*pi:pi/(f1*n_phases):cycles*pi;
subplot(2,1,1)
hold on
wave_phase = sin(x);
for iPhase = 1:cycles*n_phases

    plot(-2*f1+1+(2*f1*iPhase):(2*f1*iPhase), wave_phase(-2*f1+1+(2*f1*iPhase):(2*f1*iPhase)), 'color', c_ord(mod(iPhase,n_phases)+1,:), 'linewidth', 5)

end
xlim
axis off

%%
subplot(2,1,2)
hold on
n_phases = 10
wave_phase = sin(2.5*x);
% temp = wave_phase(1:210);
% wave_phase(1:end-210) = wave_phase(211:end);
% wave_phase(end-209:end) = temp;
for iPhase = 1:cycles*n_phases

    plot(-0.8*f1+1+(0.8*f1*iPhase):(0.8*f1*iPhase), wave_phase(-0.8*f1+1+(0.8*f1*iPhase):(0.8*f1*iPhase)), 'color', c_ord(mod(iPhase,n_phases/2.5)+1,:), 'linewidth', 5)

end
axis off