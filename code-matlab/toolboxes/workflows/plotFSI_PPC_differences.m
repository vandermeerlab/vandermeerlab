% Boilerplate code for selecting 'clean cells'

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        % do fsi stuff
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ~od.fsi_res.near_spec{iC}.flag_no_control_split
                % do fsi_stuff
                fig = figure('WindowState','maximized');
                ax = subplot(1,1,1);
                ax.YAxis.Exponent = 0;
                hold on;
                % plot the difference between Unsampled and subsampled_ppc
                keep = od.fsi_res.near_spec{iC}.mfr > 0;
                plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.ppc' - od.fsi_res.near_spec{iC}.subsampled_ppc, 'Color', 'blue')
                plot(od.fsi_res.near_spec{iC}.freqs, mean(od.fsi_res.near_spec{iC}.trialwise_ppc(keep,:)) - od.fsi_res.near_spec{iC}.subsampled_ppc, 'Color', 'green')
                plot(od.fsi_res.near_spec{iC}.freqs, mean(od.fsi_res.near_spec{iC}.trialwise_ppc(keep,:)) - od.fsi_res.near_spec{iC}.ppc', 'Color', 'magenta')
                p0 = yline(0.005, 'LineStyle', '--', 'Color','black');
                p0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                p0 = yline(-0.005, 'LineStyle', '--', 'Color','black');
                p0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                legend({'Unsampled - subsampled','Mean - subsampled', 'Mean - Unsampled'})
                clean_fsi = clean_fsi + 1;
                WriteFig(fig, strcat('FSI_PPCdif_', fsi_labels{iC}), 1);
                close;
            end
        end

        
    end
end

fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
