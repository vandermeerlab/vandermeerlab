%Script to join trial-wise and aggregated results packet
cd('D:\RandomVstrAnalysis\ft_results\');
tw_dir = 'D:\RandomVstrAnalysis\trialwise_ppc\';

rats = {'R117','R119','R131','R132'};

msn_problem = {};
fsi_problem = {};
% For now let's select the cells that have no probelmes
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        this_tw = load(strcat(tw_dir,ofiles(jdx).name));
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        %do fsi stuff
        for iC = 1:length(fsi_labels)
            % Making sure the same labels are joined
            jC = find(startsWith(this_tw.od.label(od.cell_type==2),fsi_labels{iC}));
            % Crosscheck all flags and then join
            if od.fsi_res.near_spec{iC}.flag_tooFewSpikes
                fsi_problem{length(fsi_problem)+1} = strcat(fsi_labels{iC}, ': Too Few Spikes');
            elseif od.fsi_res.near_spec{iC}.spk_count ~= this_tw.od.fsi_res.near_spec{jC}.spk_count
                fsi_problem{length(fsi_problem)+1} =strcat(fsi_labels{iC}, ': Unequal Spike Counts');
                continue;
            elseif od.fsi_res.near_spec{iC}.flag_nansts ~= this_tw.od.fsi_res.near_spec{jC}.flag_nansts
                fsi_problem{length(fsi_problem)+1} = strcat(fsi_labels{iC}, ': Nan STS not matching');
                continue;
            elseif od.fsi_res.near_spec{iC}.flag_nanppc ~= this_tw.od.fsi_res.near_spec{jC}.flag_nanppc
                fsi_problem{length(fsi_problem)+1} = strcat(fsi_labels{iC}, ': Nan PPC not matching');
                continue;
            elseif od.fsi_res.near_spec{iC}.flag_no_subsampling == this_tw.od.fsi_res.near_spec{jC}.subsampled_flag
                fsi_problem{length(fsi_problem)+1} = strcat(fsi_labels{iC}, ': Subsampling flag not matching');
                continue;
            else %join the results
                od.fsi_res.near_spec{iC}.trialwise_spk_count = this_tw.od.fsi_res.near_spec{jC}.trialwise_spk_count';
                od.fsi_res.near_spec{iC}.trialwise_ppc = this_tw.od.fsi_res.near_spec{jC}.trial_wise_ppc;
            end
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        %do msn stuff
        for iC = 1:length(msn_labels)
            % Making sure the same labels are joined
            jC = find(startsWith(this_tw.od.label(od.cell_type==1),msn_labels{iC}));
            % Crosscheck all flags and then join
            if od.msn_res.near_spec{iC}.flag_tooFewSpikes
                msn_problem{length(msn_problem)+1} = strcat(msn_labels{iC}, ': Too Few Spikes');
            elseif od.msn_res.near_spec{iC}.spk_count ~= this_tw.od.msn_res.near_spec{jC}.spk_count
                msn_problem{length(msn_problem)+1} = strcat(msn_labels{iC}, ': Unequal Spike Counts');
                continue;
            elseif od.msn_res.near_spec{iC}.flag_nansts ~= this_tw.od.msn_res.near_spec{jC}.flag_nansts
                msn_problem{length(msn_problem)+1} = strcat(msn_labels{iC}, ': Nan STS not matching');
                continue;
            elseif od.msn_res.near_spec{iC}.flag_nanppc ~= this_tw.od.msn_res.near_spec{jC}.flag_nanppc
                msn_problem{length(msn_problem)+1} = strcat(msn_labels{iC}, ': Nan PPC not matching');
                continue;
            else %join the results
                od.msn_res.near_spec{iC}.trialwise_spk_count = this_tw.od.msn_res.near_spec{jC}.trialwise_spk_count';
                od.msn_res.near_spec{iC}.trialwise_ppc = this_tw.od.msn_res.near_spec{jC}.trial_wise_ppc;
            end
        end
        fn_out = cat(2, 'D:\RandomVstrAnalysis\final_results\', ofiles(jdx).name);
        save(fn_out,'od');
    end
end
fsi_problem = fsi_problem';
msn_problem = msn_problem';
%% Save summary to output files
