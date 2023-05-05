% Script to get the PPC and STS significance

cd('D:\RandomVstrAnalysis\final_results\'); % Change this to your local machine location for results
rats = {'R117','R119','R131','R132'};
% Output is in the form of Label Sig-HFR-
sig_fsi = {};
sig_msn = {};
clean_msn = 0;
clean_fsi = 0;
min_freq = 2;  
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
                this_out = {fsi_labels{iC}, 0, 0, 0, 0};
                
                p1_p2_sts = od.fsi_res.near_p1_spec{iC}.subsampled_sts - od.fsi_res.near_p2_spec{iC}.subsampled_sts;
                p1_p2_ppc = od.fsi_res.near_p1_spec{iC}.subsampled_ppc - od.fsi_res.near_p2_spec{iC}.subsampled_ppc;
                mean_sts = mean(p1_p2_sts,1);
                sd_sts = std(p1_p2_sts);
                mean_ppc = mean(p1_p2_ppc, 1);
                sd_ppc = std(p1_p2_ppc);
                hfr_lfr_sts = od.fsi_res.near_hfr_spec{iC}.subsampled_sts - od.fsi_res.near_lfr_spec{iC}.subsampled_sts;
                hfr_lfr_ppc = od.fsi_res.near_hfr_spec{iC}.subsampled_ppc - od.fsi_res.near_lfr_spec{iC}.subsampled_ppc;
                
                % Check if sts or ppc values cross the threshold at
                x1 = find(round(od.fsi_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                foi = od.fsi_res.near_spec{iC}.freqs(x1:end);
                [d1, f1] = max(hfr_lfr_ppc(x1:end) - (mean_ppc(x1:end) + 3*sd_ppc(x1:end)));
                [d2, f2] = max((mean_ppc(x1:end) - 3*sd_ppc(x1:end)) - hfr_lfr_ppc(x1:end));
                [d3, f3] = max(hfr_lfr_sts(x1:end) - (mean_sts(x1:end) + 3*sd_sts(x1:end)));
                [d4, f4] = max((mean_sts(x1:end) - 3*sd_sts(x1:end)) - hfr_lfr_sts(x1:end));
                if d1 > 0
                    this_out{2} = foi(f1);
                end
                if d2 > 0
                    this_out{3} = foi(f2);
                end
                if d3 > 0
                    this_out{4} = foi(f3);
                end
                if d4 > 0
                    this_out{5} = foi(f4);
                end
                
                sig_fsi{clean_fsi} = this_out;
            end
        end

        %do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                %do msn_stuff
                clean_msn = clean_msn + 1;
                this_out = {msn_labels{iC}, 0, 0, 0, 0};
                
                p1_p2_sts = od.msn_res.near_p1_spec{iC}.sts - od.msn_res.near_p2_spec{iC}.sts;
                p1_p2_ppc = od.msn_res.near_p1_spec{iC}.ppc - od.msn_res.near_p2_spec{iC}.ppc;
                mean_sts = mean(p1_p2_sts,1);
                sd_sts = std(p1_p2_sts);
                mean_ppc = mean(p1_p2_ppc, 1);
                sd_ppc = std(p1_p2_ppc);
                hfr_lfr_sts = od.msn_res.near_hfr_spec{iC}.sts_vals - od.msn_res.near_lfr_spec{iC}.sts_vals;
                hfr_lfr_ppc = od.msn_res.near_hfr_spec{iC}.ppc' - od.msn_res.near_lfr_spec{iC}.ppc';
                
                % Check if sts or ppc values cross the threshold at
                x1 = find(round(od.msn_res.near_spec{iC}.freqs) > min_freq, 1, 'first');
                foi = od.msn_res.near_spec{iC}.freqs(x1:end);
                [d1, f1] = max(hfr_lfr_ppc(x1:end) - (mean_ppc(x1:end) + 3*sd_ppc(x1:end)));
                [d2, f2] = max((mean_ppc(x1:end) - 3*sd_ppc(x1:end)) - hfr_lfr_ppc(x1:end));
                [d3, f3] = max(hfr_lfr_sts(x1:end) - (mean_sts(x1:end) + 3*sd_sts(x1:end)));
                [d4, f4] = max((mean_sts(x1:end) - 3*sd_sts(x1:end)) - hfr_lfr_sts(x1:end));
                if d1 > 0
                    this_out{2} = foi(f1);
                end
                if d2 > 0
                    this_out{3} = foi(f2);
                end
                if d3 > 0
                    this_out{4} = foi(f3);
                end
                if d4 > 0
                    this_out{5} = foi(f4);
                end
                
                sig_msn{clean_msn} = this_out;
            end
        end
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
%% Shape output variables
temp = sig_fsi';
sig_fsi = cell(length(temp),5);
for i = 1:length(sig_fsi)
    for j = 1:5
        sig_fsi{i,j} = temp{i}{j};
    end
end

temp = sig_msn';
sig_msn = cell(length(temp),5);
for i = 1:length(sig_msn)
    for j = 1:5
        sig_msn{i,j} = temp{i}{j};
    end
end
