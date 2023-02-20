% Script to plot PPC by firing rate
% TODO: How to handle cases with
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'};
nbins = 8; % some fixed number
fsi_out = {};
msn_out = {};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        %do fsi stuff
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ~od.fsi_res.near_spec{iC}.flag_no_control_split
                %do fsi_stuff
                clean_fsi = clean_fsi + 1;
                nz_trials = find(od.fsi_res.near_spec{iC}.mfr > 0);
                tw_fr = od.fsi_res.near_spec{iC}.mfr(nz_trials);
                [count, edges, bin] = histcounts(tw_fr, nbins);
                fr_vals = (edges(1:end-1) + edges(2:end))/2;
                binned_ppc = zeros(length(f_list), length(count));
                for iF = 1:length(f_list)
                    f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                        round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    % take the average PPC across the frequency window of interest.
                    tw_ppc =  mean(od.fsi_res.near_spec{iC}.trialwise_ppc(nz_trials,f_idx),2);
                    for iB = 1:length(count)
                        binned_ppc(iF,iB) = mean(tw_ppc(bin == iB)); %expect nan_values
                    end
                end
                for iF = 1:length(f_list)
                    r_ppc = binned_ppc(iF,:);
                    r_fr = fr_vals;
                    n_ppc = (r_ppc - min(r_ppc))/(max(r_ppc) - min(r_ppc));
                    n_fr = (r_fr - min(r_fr))/(max(r_fr) - min(r_fr));
                    fsi_out{clean_fsi,iF} = [n_ppc;n_fr];
                end
            end
        end
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                %do msn_stuff
                clean_msn = clean_msn+1;
                nz_trials = find(od.msn_res.near_spec{iC}.mfr > 0);
                tw_fr = od.msn_res.near_spec{iC}.mfr(nz_trials);
                [count, edges, bin] = histcounts(tw_fr, nbins);
                fr_vals = (edges(1:end-1) + edges(2:end))/2;
                binned_ppc = zeros(length(f_list), length(count));
                for iF = 1:length(f_list)
                    f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                        round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    % take the average PPC across the frequency window of interest.
                    tw_ppc =  mean(od.msn_res.near_spec{iC}.trialwise_ppc(nz_trials,f_idx),2);
                    for iB = 1:length(count)
                        binned_ppc(iF,iB) = mean(tw_ppc(bin == iB)); %expect nan_values
                    end
                end
                for iF = 1:length(f_list)
                    r_ppc = binned_ppc(iF,:);
                    r_fr = fr_vals;
                    n_ppc = (r_ppc - min(r_ppc))/(max(r_ppc) - min(r_ppc));
                    n_fr = (r_fr - min(r_fr))/(max(r_fr) - min(r_fr));
                    msn_out{clean_msn,iF} = [n_ppc;n_fr];
                end
            end
        end
    end
end
%% Plot trend
fig = figure('WindowState', 'maximized');
for iF = 1:length(f_list)
    % Plot FSIs
    subplot(2,length(f_list),iF)
    hold on;
    for iC = 1:length(fsi_out)
        plot(fsi_out{iC}(1,:), 'color',c_list{iF}, 'LineWidth', 0.25);
    end
    temp_ppc = cell2mat(fsi_out(:,iF));
    plot(mean(temp_ppc(1:2:end,:),'omitnan'), 'color', 'black', 'LineWidth', 3);
    ylabel('Norm PPC')
    xlabel('Norm FR')
    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
    
    % Plot MSNs
    subplot(2,length(f_list),iF+length(f_list))
     hold on;
    for iC = 1:length(msn_out)
        plot(msn_out{iC}(1,:), 'color',c_list{iF}, 'LineWidth', 0.25);
    end
    temp_ppc = cell2mat(msn_out(:,iF));
    plot(mean(temp_ppc(1:2:end,:),'omitnan'), 'color', 'black', 'LineWidth', 3);
    ylabel('Norm PPC')
    xlabel('Norm FR')
    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
end
%% Plot heatmap
fig = figure('WindowState', 'maximized');
for iF = 1:length(f_list)
    % Plot FSIs
    subplot(2,length(f_list),iF)
    hold on;
    temp_ppc = cell2mat(fsi_out(:,iF));
    temp_ppc = temp_ppc(1:2:end,:);
    % Sort according to first ppc bin
    [~,sidx] = sort(temp_ppc(:,1));
    temp_ppc(:,:) = temp_ppc(sidx,:);
    imagesc(temp_ppc)
    ylabel('Norm PPC')
    xlabel('Norm FR')
    axis off
    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
    
    % Plot MSNs
    subplot(2,length(f_list),iF+length(f_list))
    hold on;
    temp_ppc = cell2mat(msn_out(:,iF));
    temp_ppc = temp_ppc(1:2:end,:);
    % Sort according to first ppc bin
    [~,sidx] = sort(temp_ppc(:,1));
    temp_ppc(:,:) = temp_ppc(sidx,:);
    imagesc(temp_ppc)
    ylabel('Norm PPC')
    xlabel('Norm FR')
    axis off
    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
end