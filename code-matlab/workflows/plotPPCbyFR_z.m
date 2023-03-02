% Script to plot PPC by firing rate
% TODO: First row: Raw PPC vs Raw FR, Second row: Norm PPC vs Z-Scored FR
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};

% Setting up parameters
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
c_list = {'blue', 'cyan', 'red', 'magenta', 'green'}; % make sure c_list and f_list have equal number of items
min_trial_spikes = 4; % Minimum number of spikes in a trial for it to be considered
min_trials = 10; % Minimum number of spike-count thresholded trials in a cell for it to be considered
edges = [-100,-2,-1.2,-0.4,0.4,1.2,2,100]; % The first and last bins are arbitrarily large to contain anything outside 2 standard deviations

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
                nz_trials = find(od.fsi_res.near_spec{iC}.mfr >= min_trial_spikes);
                if length(nz_trials) < min_trials
                    continue;
                end
                tw_fr = zscore(od.fsi_res.near_spec{iC}.mfr(nz_trials));
                this_mean = mean(od.fsi_res.near_spec{iC}.mfr(nz_trials));
                this_sd = std((od.fsi_res.near_spec{iC}.mfr(nz_trials)));
                [count, edges, bin] = histcounts(tw_fr, edges); 
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
                % Plot stuff
                fig = figure('WindowState', 'maximized');
                for iF = 1:length(f_list)
                    r_ppc = binned_ppc(iF,:);
                    % Getting rid of nan ppc valiues for sensible plotting
%                     keep = ~isnan(r_ppc); %getting rid of bins that had no fr values
%                     r_ppc = r_ppc(keep); 
                    % Also store normalized values
                    n_ppc = (r_ppc - min(r_ppc))/(max(r_ppc) - min(r_ppc));

                    % Plot raw values on top
                    subplot(3,length(f_list),iF)
                    plot(r_ppc, 'color', c_list{iF}, 'LineWidth', 3);
                    xticks([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
                    xticklabels({num2str(this_mean - 2*this_sd, '%.1f'), num2str(this_mean - 1.2*this_sd, '%.1f'), num2str(this_mean - 0.4*this_sd, '%.1f'), ...
                        num2str(this_mean + 0.4*this_sd, '%.1f'), num2str(this_mean + 1.2*this_sd, '%.1f'), num2str(this_mean + 2*this_sd, '%.1f')});
                    xlabel('Firing rate (Hz)')
                    ylabel('PPC')
                    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
                    
                    % Plot normalized values on bottm
                    subplot(3,length(f_list),iF+length(f_list))
                    plot(n_ppc, 'color', c_list{iF}, 'LineWidth', 3);
                    xticks([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
                    xticklabels({'-2', '-1.2', '-0.4', '0.4', '1.2', '2'})

                    xlabel('Z-Scored Firing rate')
                    ylabel('Norm PPC')
                end
                % Plot PPC and STS for comparison
                subplot(3,length(f_list), 2*length(f_list)+1:2*length(f_list)+2)
                plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_ppc, 'color', 'black')
                hold on;
                for iF = 1:length(f_list)
                    f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                        round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    area(od.fsi_res.near_spec{iC}.freqs(f_idx), od.fsi_res.near_spec{iC}.subsampled_ppc(f_idx), ...
                        'FaceColor', c_list{iF}, 'FaceAlpha', 0.5)                    
                end
                xlabel('Frequencies (Hz)');
                ylabel('PPC');
                title('ALL trials')

                subplot(3,length(f_list), 2*length(f_list)+3:2*length(f_list)+4)
                plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_sts, 'color', 'black')
                hold on;
                for iF = 1:length(f_list)
                    f_idx = find(round(od.fsi_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                        round(od.fsi_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    area(od.fsi_res.near_spec{iC}.freqs(f_idx), od.fsi_res.near_spec{iC}.subsampled_sts(f_idx), ...
                        'FaceColor', c_list{iF}, 'FaceAlpha', 0.5)                    
                end
                xlabel('Frequencies (Hz)');
                ylabel('STS');
                title('ALL trials')

                sgtitle(sprintf("FSI: %s", fsi_labels{iC}), 'Interpreter', 'None');
                WriteFig(fig,cat(2,'FSI_', fsi_labels{iC}),1);
                close;
            end
        end
        
        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                %do msn_stuff
                nz_trials = find(od.msn_res.near_spec{iC}.mfr >= min_trial_spikes);
                if length(nz_trials) < min_trials
                    continue;
                end
                tw_fr = zscore(od.msn_res.near_spec{iC}.mfr(nz_trials));
                this_mean = mean(od.msn_res.near_spec{iC}.mfr(nz_trials));
                this_sd = std((od.msn_res.near_spec{iC}.mfr(nz_trials)));
                [count, edges, bin] = histcounts(tw_fr, edges);
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
                % Plot stuff
                fig = figure('WindowState', 'maximized');
                for iF = 1:length(f_list)
                    r_ppc = binned_ppc(iF,:);
                    % Getting rid of nan ppc valiues for sensible plotting
%                     keep = ~isnan(r_ppc); %getting rid of bins that had no fr values
%                     r_ppc = r_ppc(keep); 
                    % Also store normalized values
                    n_ppc = (r_ppc - min(r_ppc))/(max(r_ppc) - min(r_ppc));

                    % Plot raw values on top
                    subplot(3,length(f_list),iF)
                    plot(r_ppc, 'color', c_list{iF}, 'LineWidth', 3);
                    xticks([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
                    xticklabels({num2str(this_mean - 2*this_sd, '%.1f'), num2str(this_mean - 1.2*this_sd, '%.1f'), num2str(this_mean - 0.4*this_sd, '%.1f'), ...
                        num2str(this_mean + 0.4*this_sd, '%.1f'), num2str(this_mean + 1.2*this_sd, '%.1f'), num2str(this_mean + 2*this_sd, '%.1f')});
                    xlabel('Firing rate (Hz)')
                    ylabel('PPC')
                    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)), 'color', c_list{iF})
                    
                    % Plot normalized values on bottm
                    subplot(3,length(f_list),iF+length(f_list))
                    plot(n_ppc, 'color', c_list{iF}, 'LineWidth', 3);
                    xticks([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
                    xticklabels({'-2', '-1.2', '-0.4', '0.4', '1.2', '2'})

                    xlabel('Z-Scored Firing rate')
                    ylabel('Norm PPC')
                end

                subplot(3,length(f_list), 2*length(f_list)+3:2*length(f_list)+4)
                plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.sts_vals, 'color', 'black')
                hold on;
                for iF = 1:length(f_list)
                    f_idx = find(round(od.msn_res.near_spec{iC}.freqs) >= f_list{iF}(1) & ...
                        round(od.msn_res.near_spec{iC}.freqs) <= f_list{iF}(2));
                    area(od.msn_res.near_spec{iC}.freqs(f_idx), od.msn_res.near_spec{iC}.sts_vals(f_idx), ...
                        'FaceColor', c_list{iF}, 'FaceAlpha', 0.5)                    
                end
                xlabel('Frequencies (Hz)');
                ylabel('STS');
                title('ALL trials')

                sgtitle(sprintf("MSN: %s", msn_labels{iC}), 'Interpreter', 'None');
                WriteFig(fig,cat(2,'MSN_', msn_labels{iC}),1);
                close;
            end
        end
    end
end
