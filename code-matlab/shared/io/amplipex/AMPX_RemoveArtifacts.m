function [ data ] = AMPX_RemoveArtifacts(data,artvalue)
%AMPX_RemoveArtifacts will take a data file in the AMPX format and change
%the artifact periods to whatever value is specified in artvalue

%% define the session type
if strcmp(data.hdr.Filename(end-7:end-4), 'post')
    session_type = 'post';
elseif strcmp(data.hdr.Filename(end-6:end-4), 'pre')
    session_type = 'pre';
else
    session_type = 'task';
end

%% find the artifacts file: should have been manually selected using FT
if exist('artifacts.mat')
    load('artifacts.mat')
else
    error('No artifacts file found in this directory.  To make one use the FieldTrip pipeline')
end


if exist('events.mat')
    load('events.mat')
else
    error('No events file found in this directory.  To make one use the FieldTrip pipeline')
end
artifact_times = artifacts.(session_type);
%% Replace the data values in during the artifact periods with the "artvalue"
%loop_num = 1; 
edge_smooth = 1*data.hdr.Fs; % one second of extra smoothing on each end of the artifact
for ii = length(artifact_times):-1:1
    for ichan = length(data.channels):-1:1;
        data.channels{1,ichan}(nearest(data.tvec, artifact_times(ii,1)/data.hdr.Fs)-edge_smooth: nearest(data.tvec, artifact_times(ii,2)/data.hdr.Fs)+edge_smooth) = artvalue;
        
    end
  %  loop_num = loop_num+1;
%     ProgressBar(loop_num/length(artifact_times), 1);
%     fprintf([num2str(floor((loop_num/length(artifact_times))*100)) '%'])
end

    for ichan = length(data.channels):-1:1;
        data.channels{1,ichan}(nearest(data.tvec, artifact_times(end,1)/data.hdr.Fs)-edge_smooth: nearest(data.tvec, data.tvec(end))) = artvalue;
        
    end
%% Since the data is a length that is not always divisable into sections the last time block is does not undergo artifact removal.  To compensate for this we need to remove the final time block which may still have artifacts that can cause problems later in the pipeline.
% if strcmp(session_type, 'task') ~= 1
    trial_length = median(diff(evt.(session_type).times(1:5)));
    for ichan = length(data.channels):-1:1;
        data.channels{1,ichan}(nearest(data.tvec, (evt.(session_type).times(end)+trial_length)/data.hdr.Fs):end) = artvalue;
        
    end
% end
artifacts.(session_type)(end+1,1) = (evt.(session_type).times(end)+trial_length);
artifacts.(session_type)(end,2) = length(data.channels{1});

%% Put the artifacts into the data structure
data.artifacts = artifacts;
data.artifacts_removed = ['Yes: Replaced with ' num2str(artvalue)];

end

