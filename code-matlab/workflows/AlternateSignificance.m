% Script to measure the alternate significance measure: i.e., are neurons
% significantly phase locked in a high firing-rate regime, compared to a
% low firing-rate regime.

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
sig2_fsi = {};
sig2_msn = {};
num_splits = 100; % This number is subject to change
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
                clean_fsi = clean_fsi+1;
                sig2_fsi{clean_fsi,1} = fsi_labels{iC};
                z_out = zeros(1, length(od.fsi_res.near_spec{iC}.freqs));
                p_out = zeros(1, length(od.fsi_res.near_spec{iC}.freqs));
                hfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc;
                lfr_ppc = od.fsi_res.near_lfr_spec{iC}.subsampled_ppc;
                p1_ppc = od.fsi_res.near_p1_spec{iC}.subsampled_ppc; 
                p2_ppc = od.fsi_res.near_p2_spec{iC}.subsampled_ppc;   
                hf_diff = hfr_ppc - lfr_ppc;
                c_diff = p1_ppc - p2_ppc;
                % Calculating the percentile for each frequency
                p_out = sum(abs(hf_diff) < abs(c_diff), 1)/num_splits; % Calculating p_val of PPC for each frequency
                % Positive sign for hfr > lfr and negative sign for hfr < lfr
                p_out(hf_diff < 0) = -1 * p_out(hf_diff < 0);
                % Calculating mean Z-score across all frequencies
                z_out = (hf_diff - mean(c_diff))./(std(c_diff));
                sig2_fsi{clean_fsi,2} = p_out;
                sig2_fsi{clean_fsi,3} = z_out;
            end
        end

        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                % do msn_stuff
                clean_msn = clean_msn+1;
                sig2_msn{clean_msn,1} = msn_labels{iC};
                z_out = zeros(1, length(od.msn_res.near_spec{iC}.freqs));
                p_out = zeros(1, length(od.msn_res.near_spec{iC}.freqs));
                hfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc';
                lfr_ppc = od.msn_res.near_lfr_spec{iC}.ppc';
                p1_ppc = od.msn_res.near_p1_spec{iC}.ppc; 
                p2_ppc = od.msn_res.near_p2_spec{iC}.ppc;   
                hf_diff = hfr_ppc - lfr_ppc;
                c_diff = p1_ppc - p2_ppc;
                % Calculating the percentile for each frequency
                p_out = sum(abs(hf_diff) < abs(c_diff), 1)/num_splits; % Calculating p_val of PPC for each frequency
                % Positive sign for hfr > lfr and negative sign for hfr < lfr
                p_out(hf_diff < 0) = -1 * p_out(hf_diff < 0);
                % Calculating mean Z-score across all frequencies
                z_out = (hf_diff - mean(c_diff))./(std(c_diff));
                sig2_msn{clean_msn,2} = p_out;
                sig2_msn{clean_msn,3} = z_out;
            end
        end
    end
end
%% Save useful variables