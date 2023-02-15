%Script to get rough count fof how many cells we have

cd('D:\RandomVstrAnalysis\ft_results\');
rats = {'R117','R119','R131','R132'};

msn_summary = table("label","too_few_spikes","spike_count","nan_in_sts","nan_in_ppc", "no_control");
fsi_summary = table("label","too_few_spikes","spike_count","nan_in_sts","nan_in_ppc", "subsampled", "no_control");
% For now let's select the cells that have no probelmes
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        %do fsi stuff
        for iC = 1:length(fsi_labels)
            this_row = {};
            if od.fsi_res.near_spec{iC}.flag_tooFewSpikes
                this_row = {fsi_labels{iC},od.fsi_res.near_spec{iC}.flag_tooFewSpikes,nan,nan,nan,nan,nan};
            elseif od.fsi_res.near_spec{iC}.flag_nanppc
                this_row = {fsi_labels{iC},od.fsi_res.near_spec{iC}.flag_tooFewSpikes,od.fsi_res.near_spec{iC}.spk_count, ...
                    od.fsi_res.near_spec{iC}.flag_nansts,od.fsi_res.near_spec{iC}.flag_nanppc, nan, ...
                    od.fsi_res.near_spec{iC}.flag_no_control_split};
            else
                this_row = {fsi_labels{iC},od.fsi_res.near_spec{iC}.flag_tooFewSpikes,od.fsi_res.near_spec{iC}.spk_count, ...
                    od.fsi_res.near_spec{iC}.flag_nansts,od.fsi_res.near_spec{iC}.flag_nanppc, ...
                    od.fsi_res.near_spec{iC}.flag_no_subsampling, od.fsi_res.near_spec{iC}.flag_no_control_split};
            end
            fsi_summary = [fsi_summary;this_row];
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
%         do msn stuff
        for iC = 1:length(msn_labels)
            this_row = {};
            if od.msn_res.near_spec{iC}.flag_tooFewSpikes
                this_row = {msn_labels{iC},od.msn_res.near_spec{iC}.flag_tooFewSpikes,nan,nan,nan,nan};
            else
                this_row = {msn_labels{iC},od.msn_res.near_spec{iC}.flag_tooFewSpikes,od.msn_res.near_spec{iC}.spk_count, ...
                    od.msn_res.near_spec{iC}.flag_nansts,od.msn_res.near_spec{iC}.flag_nanppc, ...
                    od.msn_res.near_spec{iC}.flag_no_control_split};
            end
            msn_summary = [msn_summary;this_row];
end
    end
end
%% Save summary to output files
writetable(msn_summary, 'msn_summary.xls');
writetable(fsi_summary, 'fsi_summary.xls');
