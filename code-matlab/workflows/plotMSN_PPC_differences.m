% Boilerplate code for selecting 'clean cells'

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                % do msn_stuff
                fig = figure('WindowState','maximized');
                ax = subplot(1,1,1);
                ax.YAxis.Exponent = 0;
                hold on;
                % plot the difference between Unsampled and subsampled_ppc
                keep = od.msn_res.near_spec{iC}.mfr > 0;
                plot(od.msn_res.near_spec{iC}.freqs, mean(od.msn_res.near_spec{iC}.trialwise_ppc(keep,:)) - od.msn_res.near_spec{iC}.ppc', 'Color', 'green')
                plot(od.msn_res.near_spec{iC}.freqs, 0.5*(od.msn_res.near_hfr_spec{iC}.ppc' + od.msn_res.near_lfr_spec{iC}.ppc') - od.msn_res.near_spec{iC}.ppc', 'Color', 'magenta')
                p0 = yline(0.005, 'LineStyle', '--', 'Color','black');
                p0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                p0 = yline(-0.005, 'LineStyle', '--', 'Color','black');
                p0.Annotation.LegendInformation.IconDisplayStyle = 'off';
                legend({'Mean - overall', '(HFR+LFR)/2 - overall'})
                clean_msn = clean_msn + 1;
%                 WriteFig(fig, strcat('mMSN_PPCdif_', msn_labels{iC}), 1);
                close;
            end
        end

        
    end
end

fprintf("Total number of clean mMSNs are %d.\n", clean_msn);
