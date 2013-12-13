%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPIKE SORTING WORKFLOW FOR AMPX DATA AND UMS2K %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
cd('D:\vandermeerlab\R036\R036-2013-05-29');
fname = 'R036-2013-05-29-preturn.dat';

%% load wideband data to spikesort; can be a single channel, or up to 4 channels
% using at least two is a good idea to separate noise from spikes!
channels_to_load = [10:13];
data = AMPX_loadData(fname,channels_to_load);

%% optional: get correlation matrix of decimated data to determine appropriate referencing
tic
data_ref = AMPX_loadData(fname,1:54,20); % note this will be slow -- it is loading a LOT of data!

corrMat = nan(54);

for ii = 1:54
    for jj = ii:54
   
        cc = corrcoef(data_ref.channels{ii},data_ref.channels{jj});
        corrMat(ii,jj) = cc(1,2);
        
    end
end  
  
imagesc(corrMat);
grid on;
toc

%% do the re-referencing // note this should be based on sensible channels
tic
av = AMPX_getAverage(fname,[1:8 10:15 17 19:23 26:30 32]); % also slow!
data = AMPX_reref(data,av);
toc

%% filter for spikes, and place in matrix appropriate for UMS2K
clear D;
invert_waveforms = 0; % if you have spikes bigger in the positive direction

tic
disp('filter')
for ii = length(data.channels):-1:1
    D(ii,:) = filter_for_spikes(data.channels{ii}, data.hdr.Fs);
end
if invert_waveforms, D = -D; end
toc

%% take a look
t_end = 200; % number of seconds to plot

figure(1)
clf

plot(D(:,1:t_end*data.hdr.Fs)')

%% detect spikes and align spikes
spikes = ss_default_params(data.hdr.Fs);
spikes.params.detect_method = 'manual'; % could be 'auto'
spikes.params.thresh = [-110 -110 -110 -110]; % set this carefully

tic
disp('ss_detect')
spikes = ss_detect({double(D')},spikes);
toc

tic
disp('ss_align')
spikes = ss_align(spikes);
toc

%% sort
disp('Clustering Data')

tic
disp('ss_kmeans')
spikes = ss_kmeans(spikes);

spikes = ss_energy(spikes);
spikes = ss_aggregate(spikes);
toc
%%
show_clusters(spikes,unique(spikes.assigns))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% KEY STEP: use the merge tool to assign labels to your clusters! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% further analysis steps will depend on "good unit" labels

%% merge tool
splitmerge_tool(spikes)

%% inspect (may not work for multichannel)
t_end = 50; end_idx = t_end*data.hdr.Fs;
clu = 16;
ch = 4;

subplot(211)
plot(data.tvec(1:end_idx),D(ch,1:end_idx)');
hold on;

[~,spk_idx] = histc(spikes.spiketimes(spikes.assigns == clu),data.tvec+((1./data.hdr.Fs)./2));
plot(data.tvec(spk_idx),D(ch,spk_idx)','.r');

set(gca,'XLim',[0 data.tvec(end_idx)]);

subplot(212)
plot(data.tvec(1:end_idx),data.channels{ch}(1:end_idx));
hold on;

[~,spk_idx] = histc(spikes.spiketimes(spikes.assigns == clu),data.tvec+((1./data.hdr.Fs)./2));
plot(data.tvec(spk_idx),data.channels{ch}(spk_idx)','.r');

set(gca,'XLim',[0 data.tvec(end_idx)]);

%%
compare_clusters(spikes)

%%
plot_features(spikes,1:9)

%%
outlier_tool(spikes,16)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IF SATISFIED, SAVE OUTPUT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ch_temp = sprintf('%d_',channels_to_load);

fname_out = sprintf('%s_ch%sspikes.mat',fname(1:end-4),ch_temp);
save(fname_out,'spikes');