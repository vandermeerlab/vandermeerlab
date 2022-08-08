%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PHOTOMETRY SIGNAL NON-ISOSBESTIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% From NLX (CW)
Fs = FP_data.final.Fs;
time = FP.tvec;
FP = FP.data;

%%
sessionTitle = 'M21-061R_';
time_ranges = [10, 50, 100];

for t_i = 1:length(time_ranges)
    t_range = 1:Fs*time_ranges(t_i);
    subplot(3, 1, t_i);
    plot(time(t_range), FP(t_range), 'Color', [0 0.5 0])
    title([sessionTitle, num2str(time_ranges(t_i))], 'Interpreter','none')
    ylabel('Fluorescence (dF/F)'); xlabel('Time (s)');
end

%% PETH
plot(left_end_peth.tvec, left_end_peth.data);
hold on;
plot(right_end_peth.tvec, right_end_peth.data);
legend('Left end', 'Right end')