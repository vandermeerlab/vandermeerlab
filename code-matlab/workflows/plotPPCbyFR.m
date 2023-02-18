% Script to plot PPC by firing rate
% TODO: How to handle cases with
cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
f_list = {[3 5], [7 9], [14 25], [40 65], [65 90]};
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
                nz_trials = find(od.fsi_res.near_spec{iC}.mfr > 0);
                tw_fr = od.fsi_res.near_spec{iC}.mfr(nz_trials);
                [count, edges, bin] = histcounts(tw_fr);
                fr_vals = (edges(1:end-1) + edges(2:end))/2;  
                binned_ppc = zeros(length(f_list), length(count));
                for iF = 1:length(f_list)
                    f_idx = find(od.fsi_res.near_spec{iC}.freqs >= f_list{iF}(1) & ...
                        od.fsi_res.near_spec{iC}.freqs <= f_list{iF}(2));
                    % take the average PPC across the frequency window of interest.
                    tw_ppc =  mean(od.fsi_res.near_spec{iC}.trialwise_ppc(nz_trials,f_idx),2);
                    for iB = 1:length(count)
                        binned_ppc(iF,iB) = mean(tw_ppc(bin == iB)); %expect nan_values
                    end
                end
                % Plot stuff
                fig = figure;
                for iF = 1:length(f_list)
                    r_ppc = binned_ppc(iF,:);
                    keep = ~isnan(r_ppc); %getting rid of bins that had no fr values
                    r_ppc = r_ppc(keep);
                    r_fr = fr_vals(keep);
                    % Sort ppc values and then use the sorted idx to sort the firing rates
                    [r_ppc, sidx] = sort(r_ppc);
                    r_fr = r_fr(sidx);
                    % Also store normalized values
                    n_ppc = (r_ppc - min(r_ppc))/(max(r_ppc) - min(r_ppc));
                    n_fr = r_fr/min(r_fr);

                    % Plot raw values on left
                    subplot(3,length(f_list),iF)
                    plot(r_ppc, r_fr);
                    ylabel('Firing rate (Hz)')
                    xlabel('PPC')
                    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)))
                    
                    % Plot normalized values on right
                    subplot(3,length(f_list),iF+length(f_list))
                    plot(n_ppc, n_fr);
                    ylabel('Norm Firing rate')
                    xlabel('Norm PPC')
                end
                % Plot PPC and STS for comparison
                subplot(3,length(f_list), 2*length(f_list)+1:2*length(f_list)+2)
                plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_ppc)
                xlabel('Frequencies (Hz)');
                ylabel('PPC');
                title('ALL trials')

                subplot(3,length(f_list), 2*length(f_list)+3:2*length(f_list)+4)
                plot(od.fsi_res.near_spec{iC}.freqs, od.fsi_res.near_spec{iC}.subsampled_sts)
                xlabel('Frequencies (Hz)');
                ylabel('STS');
                title('ALL trials')

                sgtitle(sprintf("FSI: %s", fsi_labels{iC}), 'Interpreter', 'None');
                WriteFig(fig,cat(2,'FSI_', fsi_labels{iC}),1);
                close;
            end
        end
        
        %do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                %do msn_stuff
                nz_trials = find(od.msn_res.near_spec{iC}.mfr > 0);
                tw_fr = od.msn_res.near_spec{iC}.mfr(nz_trials);
                [count, edges, bin] = histcounts(tw_fr);
                fr_vals = (edges(1:end-1) + edges(2:end))/2;  
                binned_ppc = zeros(length(f_list), length(count));
                for iF = 1:length(f_list)
                    f_idx = find(od.msn_res.near_spec{iC}.freqs >= f_list{iF}(1) & ...
                        od.msn_res.near_spec{iC}.freqs <= f_list{iF}(2));
                    % take the average PPC across the frequency window of interest.
                    tw_ppc =  mean(od.msn_res.near_spec{iC}.trialwise_ppc(nz_trials,f_idx),2);
                    for iB = 1:length(count)
                        binned_ppc(iF,iB) = mean(tw_ppc(bin == iB)); %expect nan_values
                    end
                end
                % Plot stuff
                fig = figure;
                for iF = 1:length(f_list)
                    r_ppc = binned_ppc(iF,:);
                    keep = ~isnan(r_ppc); %getting rid of bins that had no fr values
                    r_ppc = r_ppc(keep);
                    r_fr = fr_vals(keep);
                    % Sort ppc values and then use the sorted idx to sort the firing rates
                    [r_ppc, sidx] = sort(r_ppc);
                    r_fr = r_fr(sidx);
                    % Also store normalized values
                    n_ppc = (r_ppc - min(r_ppc))/(max(r_ppc) - min(r_ppc));
                    n_fr = r_fr/min(r_fr);

                    % Plot raw values on left
                    subplot(3,length(f_list),iF)
                    plot(r_ppc, r_fr);
                    ylabel('Firing rate (Hz)')
                    xlabel('PPC')
                    title(sprintf("%d Hz - %d Hz", f_list{iF}(1), f_list{iF}(2)))
                    
                    % Plot normalized values on right
                    subplot(3,length(f_list),iF+length(f_list))
                    plot(n_ppc, n_fr);
                    ylabel('Norm Firing rate')
                    xlabel('Norm PPC')
                end
                % Plot PPC and STS for comparison
                subplot(3,length(f_list), 2*length(f_list)+1:2*length(f_list)+2)
                plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.ppc)
                xlabel('Frequencies (Hz)');
                ylabel('PPC');
                title('ALL trials')

                subplot(3,length(f_list), 2*length(f_list)+3:2*length(f_list)+4)
                plot(od.msn_res.near_spec{iC}.freqs, od.msn_res.near_spec{iC}.sts_vals)
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
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
