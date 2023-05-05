% Script to generate a scatter plot of mean vs max PPC across various
% frequency ranges

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
% msn_twise = [];
% fsi_twise = [];

% Setting up parameters
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
min_trial_spikes = 20; % Minimum number of spikes in a trial for it to be considered
min_trials = 25; % Minimum number of spike-count thresholded trials in a cell for it to be considered
nbins = 7; % some fixed number
clean_msn = 0;
clean_fsi = 0;
fsi_mean= zeros(length(f_list),0);
fsi_max = zeros(length(f_list),0);
msn_mean= zeros(length(f_list),0);
msn_max = zeros(length(f_list),0);


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
                for iF = 1:length(f_list)
                    f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                        round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    fsi_mean(iF,clean_fsi) =  mean(od.fsi_res.near_spec{iC}.subsampled_ppc(f_idx));
                    fsi_max(iF,clean_fsi) =  max(od.fsi_res.near_spec{iC}.subsampled_ppc(f_idx));
                end
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                % do msn_stuff
                clean_msn = clean_msn + 1;
                for iF = 1:length(f_list)
                    f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                        round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    msn_mean(iF,clean_msn) =  mean(od.msn_res.near_spec{iC}.ppc(f_idx));
                    msn_max(iF,clean_msn) =  max(od.msn_res.near_spec{iC}.ppc(f_idx));
                end
            end
        end
    end
end

%%
fig = figure('WindowState', 'maximized');
for iF = 1:length(c_list)
    ax = subplot(2,length(c_list), iF);
    scatter(msn_mean(iF, :), msn_max(iF,:), 'MarkerEdgeColor', c_list{iF}, 'MarkerFaceColor', c_list{iF}, 'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha',0.5);
    xlabel('Mean PPC')
    ylabel('Max PPC')
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    title(sprintf('MSN %d Hz - %d Hz',f_list{iF}(1), f_list{iF}(2)));

    ax = subplot(2,length(c_list), iF + length(c_list));
    scatter(fsi_mean(iF, :), fsi_max(iF,:),  'MarkerEdgeColor', c_list{iF}, 'MarkerFaceColor', c_list{iF}, 'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha',0.5);
    xlabel('Mean PPC')
    ylabel('Max PPC')
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    title(sprintf('FSI %d Hz - %d Hz',f_list{iF}(1), f_list{iF}(2)));
end