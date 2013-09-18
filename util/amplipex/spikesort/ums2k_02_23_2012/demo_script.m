clear,clc
num_trials   = 20;
Fs           = 30000;
num_channels =  4;
trial_dur    = 20;

load waveform
Fr           = 20;
scaler1 = [6 3 1 4];
scaler2 = [2 6 3 1];

data = [];
for j = 1:num_trials
    data{j} = randn([Fs*trial_dur num_channels]);
end

% make waveform block
ws = [];
for j = 1:num_channels
    ws1(:,j) = w * scaler1(j);
    ws2(:,j) = w * scaler2(j);
end

% embed spikes
for j =1:num_trials
    locs = ceil( rand( [1 trial_dur * Fr] ) * (trial_dur*Fs - length(w))  );
    for k = 1:length(locs)
        if rand>.5, ws = ws1; else, ws = ws2; end
        data{j}(locs(k) + [1:length(w)], : ) = data{j}(locs(k) + [1:length(w)], : ) + ws;

    end
end

% example of filtering
Wp = [ 800  8000] * 2 / Fs;
Ws = [ 600 10000] * 2 / Fs;
[N,Wn] = buttord( Wp, Ws, 3, 20);
[B,A] = butter(N,Wn);
for j = 1:length(data)
   data2{j} = filtfilt( B, A, data{j} ); 
end

% run algorithm
spikes = ss_default_params(Fs);
spikes = ss_detect(data,spikes);
spikes = ss_align(spikes);
spikes = ss_kmeans(spikes);
spikes = ss_energy(spikes);
spikes = ss_aggregate(spikes);

% main tool
splitmerge_tool(spikes)

% stand alone outlier tool
outlier_tool(spikes)

%
% Note: In the code below, "clus", "clus1", "clus2", and "clus_list" are dummy
% variables.  The user should fill in these vlaues with cluster IDs found 
% in the SPIKES object after running the algorithm above.
%

% plots for single clusters
plot_waveforms( spikes, clus );
plot_stability( spikes, clus);
plot_residuals( spikes,clus);
plot_isi( spikes, clus );
plot_detection_criterion( spikes, clus );

% comparison plots
plot_fld( spikes,clus1,clus2);
plot_xcorr( spikes, clus1, clus2 );

% whole data plots
plot_features(spikes );
plot_aggtree(spikes);
show_clusters(spikes, [clus_list]);
compare_clusters(spikes, [clus_list]);

% outlier manipulation (see M-files for description on how to use)
spikes = remove_outliers( spikes, which ); 
spikes = reintegrate_outliers( spikes, indices, mini );

% quality metric functions
%
% Note: There are versions of these functions in the quality_measures 
% directory that have an interface that does not depend on the SPIKES
% structure.  These are for use by people who only want to use the quality
% metrics but do not want to use the rest of the package for sorting. 
% These functions have the same names as below but without the "ss_" prefix.
%
FN1 = ss_censored( spikes, clus1 );
FP1 = ss_rpv_contamination( spikes, clus1  );
FN2 = ss_undetected(spikes,clus1);
confusion_matrix = ss_gaussian_overlap( spikes, clus1, clus2 ); % call for every pair of clusters


