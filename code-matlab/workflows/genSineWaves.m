close all;
n_phases = 5;
c_ord = linspecer(n_phases);
figure
f1 = 100;
cycles = 5;
x = -cycles*pi:pi/(f1*n_phases):cycles*pi;
subplot(2,1,1)
hold on
wave_phase = sin(x);
for iPhase = 1:cycles*n_phases

    plot(-2*f1+1+(2*f1*iPhase):(2*f1*iPhase), wave_phase(-2*f1+1+(2*f1*iPhase):(2*f1*iPhase)), 'color', c_ord(mod(iPhase,n_phases)+1,:), 'linewidth', 5)

end
axis off

%%
subplot(2,1,2)
hold on
n_phases = 20;
wave_phase = sin(4*x);
temp = wave_phase(1:334);
wave_phase(1:end-334) = wave_phase(335:end);
wave_phase(end-333:end) = temp;
for iPhase = 1:cycles*n_phases

    plot(-0.5*f1+1+(0.5*f1*iPhase):(0.5*f1*iPhase), wave_phase(-0.5*f1+1+(0.5*f1*iPhase):(0.5*f1*iPhase)), 'color', c_ord(mod(mod(iPhase,n_phases)+1,5)+1,:), 'linewidth', 5)

end
axis off