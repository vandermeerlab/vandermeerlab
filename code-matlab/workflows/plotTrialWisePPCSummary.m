cd('D:\RandomVstrAnalysis\ft_trialwise_ppc');
rats = {'R117','R119','R131','R132'};

msn_ppc_sd = [];
fsi_ppc_sd = [];
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        
        fsi_labels  = od.label(od.cell_type == 2);
        %do fsi stuff
        for iC = 1:length(fsi_labels)
            if isfield(od, 'fsi_res') && isfield(od.fsi_res, 'near_spec') && ...
                   ~od.fsi_res.near_spec{iC}.flag_tooFewSpikes && ...
                   ~od.fsi_res.near_spec{iC}.flag_nanppc
               trial_count = size(od.fsi_res.near_spec{iC}.trial_wise_ppc,1);
               valid_trials = false(1, trial_count);
               for iT = 1:trial_count
                  valid_trials(iT) = ~isempty(find(od.fsi_res.near_spec{iC}.trial_wise_ppc(iT,:)));  
               end
               valid_trials = find(valid_trials);
               this_ppc = [];
               if isempty(valid_trials)
                   continue;
               end
               this_ppc = od.fsi_res.near_spec{iC}.trial_wise_ppc(valid_trials,:);
               fsi_ppc_sd = [fsi_ppc_sd; std(this_ppc)];
            end
        end
        
        msn_labels  = od.label(od.cell_type == 1);
        %do msn stuff
        for iC = 1:length(msn_labels)
            if isfield(od, 'msn_res') && isfield(od.msn_res, 'near_spec') && ...
                   ~od.msn_res.near_spec{iC}.flag_tooFewSpikes && ...
                   ~od.msn_res.near_spec{iC}.flag_nanppc
               trial_count = size(od.msn_res.near_spec{iC}.trial_wise_ppc,1);
               valid_trials = false(1, trial_count);
               for iT = 1:trial_count
                  valid_trials(iT) = ~isempty(find(od.msn_res.near_spec{iC}.trial_wise_ppc(iT,:)));  
               end
               valid_trials = find(valid_trials);
               this_ppc = [];
               if isempty(valid_trials)
                   continue;
               end
               this_ppc = od.msn_res.near_spec{iC}.trial_wise_ppc(valid_trials,:);
               msn_ppc_sd = [msn_ppc_sd; std(this_ppc)];
            end
        end
    end
end
% Plot histograms
q0 = mean(fsi_ppc_sd, 2);
q1 = mean(msn_ppc_sd, 2);
histogram(q0, 'FaceColor', 'green', 'FaceAlpha', 0.4);
hold on;
histogram(q1, 'FaceColor', 'red', 'FaceAlpha', 0.4);
