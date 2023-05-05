% Script to join trial-wise and aggregated results packet of 'clean cells'
cd('D:\RandomVstrAnalysis\final_results\');
tw_dir = 'D:\RandomVstrAnalysis\trialwise_ppc\unsampled\';
rats = {'R117','R119','R131','R132'};

% For now let's select the cells that have no problems
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        this_tw = load(strcat(tw_dir,ofiles(jdx).name));
        
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        % do fsi stuff
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ...
                    ~od.fsi_res.near_spec{iC}.flag_no_control_split
            % Making sure the same labels are joined
            jC = find(startsWith(this_tw.od.label(od.cell_type==2),fsi_labels{iC}));
            %join the results
                od.fsi_res.near_spec{iC}.trialwise_unsampled_ppc = ...
                    this_tw.od.fsi_res.near_spec{iC}.trialwise_unsampled_ppc;
            end
            fn_out = cat(2, 'D:\RandomVstrAnalysis\temp\', ofiles(jdx).name);
            save(fn_out,'od');
        end

    end
end
%% Save summary to output files
