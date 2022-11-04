cd('E:\Dropbox (Dartmouth College)\AnalysisResults\FieldTripResults\ft_trialwise_ppc');

% type = 1 for MSN, 2 for FSI
label = 'R117-2007-06-01-TT07_1'; type = 2; % Not Significant FSI
% label = 'R117-2007-06-18-TT07_2'; type = 2; % Significant FSI
% label = 'R132-2007-10-19-TT01_4'; type = 1; % Not Significant MSN
% label = 'R119-2007-07-10-TT11_3'; type = 1; % Significant MSN
fw = [5,100];

set(0,'DefaultAxesFontName','Helvetica')
set(0,'DefaultAxesFontWeight','bold')

fname = strcat(extractBefore(label,'-TT'),'_ft_spec.mat');
load(fname);
%%
if type == 2 % FSI Case
    fsi_labels  = od.label(od.cell_type == 2);
    fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
    iC = find(cellfun(@(x) strcmp(label, x), fsi_labels));
    tw_ppc = od.fsi_res.near_spec{iC}.trial_wise_ppc;
    fqs = od.fsi_res.near_spec{iC}.freqs;
else
    msn_labels  = od.label(od.cell_type == 1);
    msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
    iC = find(cellfun(@(x) strcmp(label, x), msn_labels));
    tw_ppc = od.msn_res.near_spec{iC}.trial_wise_ppc;
    fqs = od.msn_res.near_spec{iC}.freqs;
end
% get rid of trials with 0 spikes
tcount = size(tw_ppc,1);
nz_trials = false(1,tcount);
for iT = 1:tcount
    nz_trials(iT) = ~isempty(find(tw_ppc(iT,:)));
end  
tw_ppc = tw_ppc(nz_trials,:);
tw_ppc = tw_ppc(:,(fqs >= fw(1) & fqs <= fw(2)));
fqs = fqs(fqs >= fw(1) & fqs <= fw(2));

fig1  = figure('WindowState', 'maximized');
ax1 = gca(fig1);
p1 = imagesc(tw_ppc);
p2 = colorbar;
ax1.Box = 'off';
ax1.XTick = {};
ax1.TickDir = 'out';
ax1.YAxis.Label.String = 'Trials';
ax1.YAxis.FontSize = 18;
ax1.YAxis.FontWeight = 'normal';
p2.TickDirection = 'out';
p2.FontSize = 18;
p2.FontWeight = 'normal';

fig2  = figure('WindowState', 'maximized');
ax2 = gca(fig2);
p3 = plot(fqs, std(tw_ppc));
p3.LineWidth = 5;
ax2.Box = 'off';
ax2.TickDir = 'out';
ax2.YAxis.Label.String = 'SD';
ax2.YAxis.FontSize = 18;
ax2.YAxis.FontWeight = 'normal';
ax2.XAxis.Label.String = 'Frequency (Hz)';
ax2.XAxis.FontSize = 18;
ax2.XAxis.FontWeight = 'normal';
ax2.YAxis.Exponent = 0;
ax2.XLim = fw;
ax2.XTick = [5 10 20 30 40 50 60 70 80 90 100];
ax2.YLim = [-0.1 0.5];
% Put breakpoint below for this to work
% ax2.Position(3) = ax1.Position(3); % Hacky way to align plots
 %%
 function ndata = normdata(data)
    if (sum(isnan(data)) == length(data))
        ndata = data;
    else
        maxd = max(data);
        mind = min(data);
        ndata = (data - mind)/(maxd-mind);
    end
end
