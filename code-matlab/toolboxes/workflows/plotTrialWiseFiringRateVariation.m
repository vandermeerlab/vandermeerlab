% Plotting the distribution of variance in trial-wise normalized firing-rates

cd('D:\RandomVstrAnalysis\tria\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
fsi_fr_sd = [];
msn_fr_sd = [];
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
                clean_fsi = clean_fsi + 1;
                tw_mfr = od.fsi_res.near_spec{iC}.mfr;
                tw_mfr = tw_mfr(tw_mfr > 0); % Getting rid of 0 spike trials
                tw_mfr = (tw_mfr - min(tw_mfr))/(max(tw_mfr) - min(tw_mfr)); % Normalize values between 0-1
                fsi_fr_sd = [fsi_fr_sd, std(tw_mfr)];
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                % do msn_stuff
                clean_msn = clean_msn + 1;
                tw_mfr = od.msn_res.near_spec{iC}.mfr;
                tw_mfr = tw_mfr(tw_mfr > 0); % Getting rid of 0 spike trials
                tw_mfr = (tw_mfr - min(tw_mfr))/(max(tw_mfr) - min(tw_mfr)); % Normalize values between 0-1
                msn_fr_sd = [msn_fr_sd, std(tw_mfr)];
            end
        end
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);

%% Plot histograms
fig = figure('WindowState', 'maximized');
ax1 = subplot(2,1,1);
h1 = histogram(fsi_fr_sd, 0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(msn_fr_sd,0:0.005:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax1.XAxis.FontSize = 20;
ax1.YAxis.FontSize = 20;
ax1.TickDir = 'out';
ax1.YTick = [0 0.5 1];
ax1.XTick = [0:0.1:0.5];
ax1.XLabel.String = 'Norm Firing Rate SD';
ax1.YLabel.String = 'Proportion';
ax1.XLim = [0 0.5];
ax1.YLim = [0 1];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 20;
leg.FontWeight = 'bold';
ax1.Box = 'off';
