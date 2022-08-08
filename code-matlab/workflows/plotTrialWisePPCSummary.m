cd('D:\RandomVstrAnalysis\ft_trialwise_ppc');
% cd('/Users/manishm/Dropbox (Dartmouth College)/AnalysisResults/FieldTripResults/ft_trialwise_ppc');
rats = {'R117','R119','R131','R132'};

msn_ppc_sd = [];
fsi_ppc_sd = [];
missing_msn = [];

load('./cellsOfInterest.mat');
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
            if sum(contains(fsi_counter, fsi_labels(iC))) == 1
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
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
%         do msn stuff
        for iC = 1:length(msn_labels)
            if sum(contains(msn_counter, msn_labels(iC))) == 1
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
%% Plot histograms
close all;
q0 = mean(fsi_ppc_sd, 2);
q1 = mean(msn_ppc_sd, 2);
fig = figure('WindowState', 'maximized');
h1 = histogram(q0,0:0.05:1, 'FaceColor', 'green', 'FaceAlpha', 0.4);
hold on;
h2 = histogram(q1,0:0.05:1, 'FaceColor', 'red', 'FaceAlpha', 0.4);
ax = gca(fig);
leg = legend({'FSI', 'MSN'});
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XLabel.String = 'Mean of SDs over all frequencies';
ax.XLabel.FontSize = 24;
leg = legend({'MSN', 'FSI'});
leg.FontName = 'Arial';
leg.FontSize = 20;
leg.FontWeight = 'bold';
box off;
%%
title('Distribution of standard deviations of ppc across trials, averaged over all the frequencies')

%% Helper functions
% function to normalize data
function ndata = normdata(data)
    if (sum(isnan(data)) == length(data))
        ndata = data;
    else
        maxd = max(data);
        mind = min(data);
        ndata = (data - mind)/(maxd-mind);
    end
end
