% Script to join trial-wise and aggregated results packet of 'clean cells'
cd('D:\RandomVstrAnalysis\ft_results\');
tw_dir1 = 'D:\RandomVstrAnalysis\trialwise_ppc\';
tw_dir2 = 'D:\RandomVstrAnalysis\trialwise_ppc\noThresh\';
rats = {'R117','R119','R131','R132'};

% For now let's select the cells that have no probelmes
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        this_tw1 = load(strcat(tw_dir1,ofiles(jdx).name));
        this_tw2 = load(strcat(tw_dir2,ofiles(jdx).name));
        
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        % do fsi stuff
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ...
                    ~od.fsi_res.near_spec{iC}.flag_no_control_split
            % Making sure the same labels are joined
            jC = find(startsWith(this_tw1.od.label(od.cell_type==2),fsi_labels{iC}));
            %join the results
                od.fsi_res.near_spec{iC}.thresh_trialwise_spk_count = this_tw1.od.fsi_res.near_spec{jC}.trialwise_spk_count';
                od.fsi_res.near_spec{iC}.thresh_trialwise_ppc = this_tw1.od.fsi_res.near_spec{jC}.trialwise_ppc;
                od.fsi_res.near_spec{iC}.no_thresh_trialwise_spk_count = this_tw2.od.fsi_res.near_spec{jC}.trialwise_spk_count';
                od.fsi_res.near_spec{iC}.no_thresh_trialwise_ppc = this_tw2.od.fsi_res.near_spec{jC}.trialwise_ppc;
            end
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        %do msn stuff
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ...
                    ~od.msn_res.near_spec{iC}.flag_no_control_split
                % Making sure the same labels are joined
                jC = find(startsWith(this_tw1.od.label(od.cell_type==1),msn_labels{iC}));
                %join the results
                od.msn_res.near_spec{iC}.thresh_trialwise_spk_count = this_tw1.od.msn_res.near_spec{jC}.trialwise_spk_count';
                od.msn_res.near_spec{iC}.thresh_trialwise_ppc = this_tw1.od.msn_res.near_spec{jC}.trialwise_ppc;
                od.msn_res.near_spec{iC}.no_thresh_trialwise_spk_count = this_tw2.od.msn_res.near_spec{jC}.trialwise_spk_count';
                od.msn_res.near_spec{iC}.no_thresh_trialwise_ppc = this_tw2.od.msn_res.near_spec{jC}.trialwise_ppc;
            end
            fn_out = cat(2, 'D:\RandomVstrAnalysis\final_results\', ofiles(jdx).name);
            save(fn_out,'od');
        end
    end
end
%% Save summary to output files
