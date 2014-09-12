function [fnBase, nFiles] = WriteKKwikFeatureFile(fnBase, features, varargin)

% [fnBase, nFiles] = WriteKKwikFeatureFile(fnBase, features, varargin)
%
% Parameters 
%    maxSpikesBeforeSplit; 0 means don't split
%    spikes; [] means use all spikes

MCD = MClust.GetData();

% -------------------
% Parms

maxSpikesBeforeSplit = 0;  % use 0 for don't split
spikes = [];
FeatureTimestamps = MCD.FeatureTimestamps;

process_varargin(varargin);


% -------------------

nF = length(features);
if isempty(spikes)
    nSpikes = length(FeatureTimestamps);
else
    nSpikes = length(spikes);
end

FD = nan(nF, nSpikes);
% count features
for iF = 1:nF
    if isempty(spikes)
        FD(iF,:) = features{iF}.GetData();
    else
        FD0 = features{iF}.GetData();
        FD(iF,:) = FD0(spikes);
    end
    if ~features{iF}.normalizedYN % KKwik needs it normalized
        FD(iF,:) = (FD(iF,:)-mean(FD(iF,:)))/std(FD(iF,:));
    end
end

if ~isempty(maxSpikesBeforeSplit) && maxSpikesBeforeSplit>0 && nSpikes > maxSpikesBeforeSplit
    % split into files
    spikeSplit = 1:maxSpikesBeforeSplit:(nSpikes) ;
    nFiles = length(spikeSplit);
    for iFile = 1:(nFiles-1)
        
        fn = [fnBase '.fet.' num2str(iFile)];
        dlmwrite(fn, sprintf('%d', nF));
        dlmwrite(fn, FD(:,spikeSplit(iFile):(spikeSplit(iFile+1)-1))', '-append', 'delimiter', ' ');
        
    end
else
    iFile = 1;
    
    fn = [fnBase '.fet.' num2str(iFile)];
    dlmwrite(fn, sprintf('%d', nF));
    dlmwrite(fn, FD', '-append', 'delimiter', ' ');
    
end

nFiles = max(iFile);

end % WriteKKwikFeatures
