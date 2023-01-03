cd('E:\Dropbox (Dartmouth College)\AnalysisResults\FieldTripResults\ft_trialwise_ppc');
rats = {'R117','R119','R131','R132'};

msn_ppc_sd = [];
fsi_ppc_sd = [];
missing_msn = [];

set(0,'DefaultAxesFontName','Helvetica');
set(0,'DefaultAxesFontWeight','bold');

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
h1 = histogram(q0,0:0.05:1, 'Normalization', 'probability', 'FaceColor', 'green', 'FaceAlpha', 1);
hold on;
h2 = histogram(q1,0:0.05:1, 'Normalization', 'probability', 'FaceColor', 'red', 'FaceAlpha', 0.6);
ax = gca(fig);
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
ax.TickDir = 'out';
ax.YTick = [0 0.5 1];
ax.XTick = [0 0.3 0.5];
ax.XLabel.String = 'Mean of SDs over all frequencies';
ax.YLabel.String = 'Proportion';
% ax.XLabel.FontSize = 40;
ax.XLim = [0, 0.5];
leg = legend({'FSI', 'MSN'});
leg.FontName = 'Helvetica';
leg.FontSize = 40;
leg.FontWeight = 'bold';
box off;
% title('Distribution of Average SD', 'FontSize', 40)

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
