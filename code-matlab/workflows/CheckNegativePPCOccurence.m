% Code to check how often does negative PPC values occur as a function ofnumber of spikes

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
msn_sc = [];
msn_nf = [];
fsi_nf = [];
fsi_sc = [];

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
                fsi_sc = [fsi_sc, od.fsi_res.near_spec{iC}.spk_count];
                fsi_nf = [fsi_nf, ...
                    sum(od.fsi_res.near_spec{iC}.subsampled_ppc < 0)/length(od.fsi_res.near_spec{iC}.subsampled_ppc)];
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                % do msn_stuff
                clean_msn = clean_msn + 1;
                msn_sc = [msn_sc, od.msn_res.near_spec{iC}.spk_count];
                msn_nf = [msn_nf, ...
                    sum(od.msn_res.near_spec{iC}.ppc < 0)/length(od.msn_res.near_spec{iC}.ppc)];
            end
        end
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
%%
fig = figure;
subplot(2,1,1);
scatter(msn_sc, msn_nf, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', ...
    'MarkerFaceAlpha', 0.5, 'SizeData',50);
xlabel('Spike Count');
ylabel('Proportion of negative PPC')
title('MSN')

subplot(2,1,2)
scatter(fsi_sc, fsi_nf, 'MarkerFaceColor', 'blue', 'MarkerFaceAlpha', 0.5, 'SizeData',50);
xlabel('Spike Count');
ylabel('Proportion of negative PPC')
title('FSI')
