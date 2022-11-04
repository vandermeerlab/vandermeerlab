cd('E:\Dropbox (Dartmouth College)\AnalysisResults\FieldTripResults\ft_results');
rats = {'R117','R119','R131','R132'};

msn_labels = {};
fsi_labels = {};
freqs = [];
msn_ppc = [];
fsi_ppc= [];
freqs = [];
freq_flag = false;
load('cellsOfInterest.mat');

% Setting default plot parameters
set(0,'DefaultAxesFontName','Helvetica');
set(0,'DefaultAxesFontWeight','bold');
%%
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        
        % do fsi stuff
        this_labels  = od.label(od.cell_type == 2);
        this_labels = cellfun(@(x) extractBefore(x, '.t'), this_labels, 'UniformOutput', false);
        for iC = 1:length(this_labels)
            if sum(contains(fsi_counter, this_labels(iC))) == 1
                if ~freq_flag
                   freqs = od.fsi_res.near_spec{iC}.freqs;
                   freq_flag = true;
                end
                fsi_labels{length(fsi_labels)+1} = this_labels{iC};
                fsi_ppc = [fsi_ppc; od.fsi_res.near_spec{iC}.subsampled_ppc];
            end
        end
        
        % do msn stuff
        this_labels  = od.label(od.cell_type == 1);
        this_labels = cellfun(@(x) extractBefore(x, '.t'), this_labels, 'UniformOutput', false);
        for iC = 1:length(this_labels)
            if sum(contains(msn_counter, this_labels(iC))) == 1
                if ~freq_flag
                   freqs = od.msn_res.near_spec{iC}.freqs;
                   freq_flag = true;
                end
                msn_labels{length(msn_labels)+1} = this_labels{iC};
                msn_ppc = [msn_ppc; od.msn_res.near_spec{iC}.ppc'];
            end
        end
    end      
end

%% Create MSN Histograms for various frequency bands of interest
f_list = {[2 5], [6 10], [25 55], [65 100]};
axes = cell(size(f_list));
fig = figure('WindowState', 'Maximized');
for iF = 1:length(f_list)
    axes{iF} = subplot(2, 2, iF);
    f_start = f_list{iF}(1);
    f_stop = f_list{iF}(2);
    subset_freqs = (freqs >= f_start) & (freqs <= f_stop);
%     msn_subset = max(msn_ppc(:,subset_freqs), [], 2);  
    msn_subset = mean(msn_ppc(:,subset_freqs), 2);  
    h1 = histogram(msn_subset, [0:0.00025:0.05,1],'Normalization', 'probability', 'FaceColor', 'black', 'FaceAlpha', 1); 
    h1.EdgeAlpha = 1;
    h1.EdgeColor = [1 1 1];
    title(sprintf('%d Hz - %d Hz', f_start, f_stop), 'FontSize', 25);
    grid on;
    axes{iF}.XLim = [0 0.051];
    axes{iF}.YLim = [0 0.1];
    axes{iF}.XAxis.FontSize = 18;
    axes{iF}.XAxis.FontWeight = 'normal';
    axes{iF}.YAxis.FontSize = 18;
    axes{iF}.YAxis.FontWeight = 'normal';
    axes{iF}.TickDir = 'out';
    axes{iF}.XAxis.Label.String = 'PPC';
    axes{iF}.YAxis.Label.String = 'Proportion';
    axes{iF}.Box = 'off';
end
%Put breakpoint and then run
for iF = 1:length(f_list)
    axes{iF}.XAxis.TickLabels{end} = '>0.05';
    axes{iF}.YAxis.TickLabels{end} = '>0.1';
end

%%
f_list = {[2 5], [6 10], [25 55], [65 100]};
fig = figure('WindowState', 'Maximized');
for iF = 1:length(f_list)
    axes{iF} = subplot(2,2, iF);
    f_start = f_list{iF}(1);
    f_stop = f_list{iF}(2);
    subset_freqs = (freqs >= f_start) & (freqs <= f_stop);
%     fsi_subset = max(fsi_ppc(:,subset_freqs), [], 2);  
    fsi_subset = mean(fsi_ppc(:,subset_freqs), 2); 
    h1 = histogram(fsi_subset, [0:0.00025:0.05,1], 'Normalization', 'probability', 'FaceColor', 'black', 'FaceAlpha', 1); 
    h1.EdgeAlpha = 1;
    h1.EdgeColor = [1 1 1];
    title(sprintf('%d Hz - %d Hz', f_start, f_stop), 'FontSize', 25);
    grid on;
    axes{iF}.XLim = [0 0.051];
    axes{iF}.YLim = [0 0.1];
    axes{iF}.XAxis.FontSize = 18;
    axes{iF}.XAxis.FontWeight = 'normal';
    axes{iF}.YAxis.FontSize = 18;
    axes{iF}.YAxis.FontWeight = 'normal';
    axes{iF}.TickDir = 'out';
    axes{iF}.XAxis.Label.String = 'PPC';
    axes{iF}.YAxis.Label.String = 'Proportion';
    axes{iF}.Box = 'off';
end
%Put breakpoint and then run
for iF = 1:length(f_list)
    axes{iF}.XAxis.TickLabels{end} = '>0.05';
    axes{iF}.YAxis.TickLabels{end} = '>0.1';
end

%% testing pctrile function and plotting
figure;
req_msn_label = 'R117-2007-06-09-TT03_4';
label_idx = find(strcmp(msn_labels, req_msn_label));
test = get_ptile(msn_ppc, label_idx);
bar(1:99/91:99.6, test)


%% Helper functions
% function to return frequency wise percentile values for input cells
function output = get_ptile(data, indices)
    [nrows, ncols] = size(data);
    output = zeros(length(indices), ncols);
    for i = 1:length(indices)
        idx = indices(i);
        for iR = 1:ncols
            output(i,iR) = (length(find(data(:,iR)< data(idx,iR))))/nrows * 100;
        end         
    end
end
