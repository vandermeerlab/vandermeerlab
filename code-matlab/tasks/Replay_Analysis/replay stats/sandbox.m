%% Get data
originalFolder = pwd;
cfg = [];
% cfg.rats = {'R042','R044','R050'};
cfg.rats = {'R050'};
% cfg.restrict = 'pre';
cfg.NAU = 1;

% profile on
data = GenSWRstats_all(cfg);
simdata = SimSWRstats_all(cfg,data);
% profile viewer

cd(originalFolder)
%% Unpack stats across sessions

% IFstats.FR_cell = [];
% IFstats.FR_field = [];
% IFstats.ISI = [];
% sess = 1;
% session.avgspk{sess} = [];

data_in = [data simdata];
stats(2).tlen = [];
stats(2).mean = [];
stats(2).var = [];
stats(2).prop = [];
stats(2).avgspk = [];
stats(2).ISI = [];
stats(2).pathlen = [];
    
for iD=1:2
    for iR=1:length(cfg.rats)
        for iT=1:2
            rn = cfg.rats{iR};
            for iFD=1:length(data_in(iD).R(iT).(rn))
                curr_R = data_in(iD).R(iT).(rn){1,iFD};

                [nCells,nSWR] = size(curr_R);
                if isempty(curr_R)
                    continue 
                end

                % get mean and variance across SWRs (for each unit)
                stats(iD).mean = vertcat(stats(iD).mean,mean(curr_R,2));
                stats(iD).var = vertcat(stats(iD).var,var(curr_R,0,2));

                % get length of each SWR event
                stats(iD).tlen = vertcat(stats(iD).tlen,data_in(iD).SWR_size(iT).(rn){iFD});

                % proportion of SWR each cell is active
                stats(iD).prop = vertcat(stats(iD).prop,sum(curr_R>0,2)/nSWR*100.0);

                % mean spike count of place cells across SWR events
                curr_R(curr_R==0) = nan; %SWR events with no spikes are omitted from the average
                stats(iD).avgspk = vertcat(stats(iD).avgspk,nanmean(curr_R,2));
%                 session.avgspk{sess} = nanmean(curr_R,2);
%                 sess = sess + 1;

                % get ISIs
                stats(iD).ISI = vertcat(stats(iD).ISI,data_in(iD).ISI_SWR(iT).(rn){iFD});
%                 IFstats.ISI = vertcat(IFstats.ISI,data_in.ISI_IF(iT).(rn){iFD});
% 
%                 % get infield firing rates
%                 IFstats.FR_cell = vertcat(IFstats.FR_cell,data_in.IFFR(iT).(rn){iFD});
%                 IFstats.FR_field = vertcat(IFstats.FR_field,data_in.IFFR_field(iT).(rn){iFD});

                % get path lengths
                stats(iD).pathlen = vertcat(stats(iD).pathlen,data_in(iD).path_len(iT).(rn){iFD});
                
            end %iterate session
            
        end %iterate left vs right
        
    end %iterate rats
    
end %iterate empirical vs simulation


%% Plot SWR data vs simulated data for varying statistics

data_in = stats;
data_name = {'Empirical','Model'};

% participation ratio
figure('units','normalized','outerposition',[0,0,1,1]); 
for i=1:length(data_in)
    subplot(1,2,i)
    [n,x] = hist(data_in(i).prop,20);
    bar(x,n,'hist')
    set(gca,'FontSize',12);
    xlabel('% SPWR participation')
    ylabel('Number of cells')
    title(data_name(i));
end
print(gcf,'-dpng','-r300','Cell participation comparison');
clear x n i


% SWR spike count distributions for place cells
figure('units','normalized','outerposition',[0,0,1,1]);
for i=1:length(data_in)
    subplot(1,2,i)
    [n,x] = hist(data_in(i).avgspk);
    bar(x,n./sum(n)*100,'hist')
    set(gca,'FontSize',12); 
    xlabel('Mean spike count')
    ylabel('% of place fields')
    title(data_name(i));
end
print(gcf,'-dpng','-r300','Spike count distribution comparison')
clear x n


% ISI distributions
dt = 0.005; % in s, because spike times are in s
tmax = 0.2;
isi_edges = 0:dt:tmax; % bin edges for ISI histogram
isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
figure('units','normalized','outerposition',[0,0,1,1]); 
for i=1:length(data_in)
    subplot(1,2,i)
    isih = histc(data_in(i).ISI,isi_edges);
    isih = isih/sum(isih);

    bar(isi_centers*1000,isih(1:end-1)); % remember to ignore the last bin of histc() output
    set(gca,'FontSize',12); 
    xlabel('ISI (ms)'); ylabel('percent of isi'); grid on;
    title(data_name(i))
end
print(gcf,'-dpng','-r300','ISI distribution comparison')
clear dt histdata isi_centers isi_edges isih tmax


% Get SWR interval distribution
dt = 0.01; % in s, because spike times are in s
figure('units','normalized','outerposition',[0,0,1,1]); 
for i=1:length(data_in)
    tmax = max(data_in(i).tlen);
    isi_edges = 0:dt:tmax; % bin edges for ISI histogram
    isi_centers = isi_edges(1:end-1)+dt/2; % for plotting
    isih = histc(data_in(i).tlen,isi_edges);
    isih = isih/sum(isih);

    subplot(1,2,i)
    bar(isi_centers*1000,isih(1:end-1));
    set(gca,'FontSize',12); xlabel('SWR interval (ms)'); ylabel('percent intervals'); grid on;
    title(data_name(i))
end
print(gcf,'-dpng','-r300','Event length distribution comparison')
clear dt i tmax isi_edges isi_centers isih

% Fano Factor comparison
figure('units','normalized','outerposition',[0,0,1,1]); 
for i=1:length(data_in)
    poiss_mean = [];
    poiss_var = [];
    for iC=1:length(data_in(i).mean)
        lambda = data_in(i).mean(iC);
        rand_sample = poissrnd(lambda,1,1000);
        poiss_mean = vertcat(poiss_mean,mean(rand_sample,2));
        poiss_var = vertcat(poiss_var,var(rand_sample,0,2));
    end

%     % plot ratio as scatter plot
%     [p,s] = polyfit(data_in(i).mean,data_in(i).var,1);
%     [y,d] = polyval(p,data_in(i).mean,s);
    subplot(1,2,i)
    scatter(data_in(i).mean,data_in(i).var,'x'); hold on;
    scatter(poiss_mean,poiss_var,'xr')
    edgeval = [floor(min(data_in(i).var)) ceil(max(data_in(i).var))];
    plot([edgeval(1),edgeval(2)],[edgeval(1),edgeval(2)],'k:')
    set(gca,'FontSize',12); 
    xlabel('mean'); ylabel('variance')
    title(data_name(i))
%     clear p s y d iC
end
legend({'data','poisson generated','ideal ratio'},'FontSize',8)
print(gcf,'-dpng','-r300','fano factor distribution comparison')
clear i iC

% path length comparison
figure('units','normalized','outerposition',[0,0,1,1]); 
for i=1:length(data_in)
    subplot(1,2,i)
    [n,x] = hist(data_in(i).pathlen,20);
    bar(x,n,'hist')
    set(gca,'FontSize',12);
    xlabel('path length')
    ylabel('Number of cells')
    title(data_name(i));
end
print(gcf,'-dpng','-r300','Replay path length comparison');
clear x n i


%% autocorrelation

test = restrict2(S_pc(1),evt);
pair = combnk(1:length(test.t),2);
xcorr_all = zeros(length(pair),1001);
for iP=1:length(pair)
    [xcorr,xbin] = ccf([],test.t{pair(iP,1)},test.t{pair(iP,2)});
    xcorr_all(iP,:) = xcorr;
end

bar(xbin,mean(xcorr_all))
xlabel('Lag (s)')
ylabel('mean xcorr')
% print(gcf,'-dpng','-r300','cand_xcorrs');

%% Get correlations
x = [];
y1 = []; y2=y1;
for iT=1:2
    for iR=1:3
        rn = rats{iR};
        x = vertcat(x,IFFR(iT).(rn));
        y1 = vertcat(y1,SWR_prop(iT).(rn));
        y2 = vertcat(y2,SWR_avgspk(iT).(rn));
    end
end

[rho1,pval1] = corr(x,y1)
[rho2,pval2] = corr(x,y2)

%% Plot firing rate distributions for place cells in their fields
figure('units','normalized','outerposition',[0,0,1,1]); 
[n,x] = hist(IFstats.FR_field);
bar(x,n./sum(n)*100,'hist')
xlabel('Mean in-field firing rate (Hz)')
ylabel('% of place fields')
% print(gcf,'-dpng','-r300','Firing rate distributions for place fields')
clear x n

%% Plot in-field firing rate vs SWR participation
figure('units','normalized','outerposition',[0,0,1,1]); 
hold on;
col = 'rb';
shape = 'o+d';
for iR=1:length(rats)
    rn = rats{iR};
    for iT=1:2
        scatter(data_in.IFFR(iT).(rn),SWR_prop(iT).(rn),20,col(iT),shape(iR));
    end
end

suptitle('SWR participation vs in-field firing rate')
ylabel('SWR Particpation')
xlabel('In-field Firing Rate')
legend({'R042 (vs left)','vs right','R044 (vs left)','vs right','R050 (vs left)','vs right'},'Location','NorthWest')
% print(gcf,'-dpng','-r300','IFFR_SWRp_scatter');

%% plot infield FR vs SWR participation for different sessions
nsp = numSubplots(32);
iP = 1;
figure;
for iT=1:2
    for iR=1:3
        rn = cfg.rats{iR};
        for iFD=1:length(data_in.R(iT).(rn))
            subplot(nsp(1),nsp(2),iP)
            curr_R = data_in.R(iT).(rn){1,iFD};
            if isempty(curr_R)
                continue
            end
            
            temp_FR = data_in.IFFR(iT).(rn){1,iFD};
            temp_prop = sum(curr_R>0,2)/length(curr_R)*100.0;
%             [rho,pval] = corr(temp_FR,temp_prop);
            scatter(temp_FR,temp_prop);
%             legend(sprintf('%0.4f',rho))
            iP = iP + 1;
        end
    end
end
clear nsp iT iR rn iFD curr_R temp_FR temp_prop iP



%% Plot in-field firing rate vs SWR spike count

figure('units','normalized','outerposition',[0,0,1,1]); 
hold on;
col = 'rb';
shape = 'o+d';
for iR=1:length(rats)
    rn = rats{iR};
    for iT=1:2
        scatter(IFFR(iT).(rn),SWR_avgspk(iT).(rn),20,col(iT),shape(iR));
    end
end
suptitle('mean SWR spike count vs in-field firing rate')
ylabel('SWR spike count')
xlabel('In-field firing rate')
legend({'R042 (vs left)','vs right','R044 (vs left)','vs right','R050 (vs left)','vs right'},'Location','NorthWest')
print(gcf,'-dpng','-r300','IFFR_SWRspkcount_scatter');

%% Get spike count histograms across SWRs

curr_R = data_in.R(2).R050{1,1};
curr_R(curr_R==0) = nan;
[nCells,nSWR] = size(curr_R);
[nsp,~] = numSubplots(nCells);
dx = 1; % in s, because spike times are in s
xmax = 10;
hedges = 0:dx:xmax; % bin edges for ISI histogram
hcenters = hedges(1:end-1)+dx/2; % for plotting

figure('units','normalized','outerposition',[0,0,1,1]); 
for iC=1:nCells
    subplot(nsp(1),nsp(2),iC)
    curr_cell = curr_R(iC,:);
    n = hist(curr_cell,hedges);
    bar(hcenters,n(1:end-1));
    set(gca,'XLim',[0,10],'XTick',0:10)
end
clear curr_R nCells nSWR nsp dx xmax hedges hcenters iC curr_cell n



%%










