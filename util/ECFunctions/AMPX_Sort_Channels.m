function [ data ] = AMPX_Sort_Channels(data, varargin )
%AMPX_Sort_Channel:
% Input a data det from the AMPX_loadData and the type of probe used
% (either 'Buz64' or 'A8x8').  This will then reorder the channels based on
% the layout of the probe.  This can also be modified to cope with
% recording session that used the reverse connection layout between the
% headstage and the probe board.

% INPUTS:
%
% data: [1x1 struct] - this is the output of AMPX_loadData containing the
% hdr, channels, tvec, and labels.
% Probe: [string] this is the name of the probe type used (either 'Buz64'
% or 'A8x8'
% Connection: [1 x 1 logical] if the probe was connectect
% 'Omnetics-Omnetics' (standard) then this 'Connection' = 1.  If it is
% reversed then it is 0.
%
%
% OUTPUTS:
%
% data: [1 x 1 struct] with the channel labels corresponding to the correct
% location.

%% initialize the variables:

probe = 'Buz64';
Connection = 1;  % this is the standard Omnetics-Omnetics connection

extract_varargin;

%% Probe impedance maps
if strcmp(data.hdr.Filename(5:8),'R036')==1
    data.Imp_channels = [.356 .316 .373 .292 .330 .384 .276 .139 .284 .285 .336 2.338 .3 .355 .353 .348 .292 .386 .324 .302 .322 .284 .294 .375 .319 .301 .376 .329 .346 .323 .355 .305 .349 .353 .331 .236 .364 .383 .421 .108 .343 .317 .279 .339 2.592 2.594 2.652 2.630 .312 .386 .323 .252 .360 .351 .310 .320 .495 .258 .280 .305 2.671 2.645 2.607 2.767];
    data.bad_channels = [1:8,33,34];
    
end

%%  Create a probe conversion script for the Buzsaki 64 probe
if strcmp(probe,'Buz64')
    Input_labels = [33 40 36 39 35 38 36 37 41 48 42 47 43 46 44 45 49 56 50 55 51 54 52 53 57 64 58 63 59 62 60 61 1 8 2 7 3 6 4 5 9 16 10 15 11 14 12 13 17 24 18 23 17 22 20 21 25 32 26 31 27 30 28 29];
    Loc_labels = [1 8 2 7 3 6 4 5 9 16 10 15 11 14 12 13 17 24 18 23 17 22 20 21 25 32 26 31 27 30 28 29 33 40 36 39 35 38 36 37 41 48 42 47 43 46 44 45 49 56 50 55 51 54 52 53 57 64 58 63 59 62 60 61];
    
    %     if Connection == 1;
    %         input_labels = [1 8 2 7 3 6 4 5 9 16 10 15 11 14 12 13 17 24 18 23 17 22 20 21 25 32 26 31 27 30 28 29 33 40 36 39 35 38 36 37 41 48 42 47 43 46 44 45 49 56 50 55 51 54 52 53 57 64 58 63 59 62 60 61];
    %         Loc_labels = [1 8 2 7 3 6 4 5 9 16 10 15 11 14 12 13 17 24 18 23 17 22 20 21 25 32 26 31 27 30 28 29 33 40 36 39 35 38 36 37 41 48 42 47 43 46 44 45 49 56 50 55 51 54 52 53 57 64 58 63 59 62 60 61];
    %     end
    
    data.labels = Loc_labels;
    temp_channels = cell(1,64);
    for iChan = 1:size(data.channels,2)
        temp_channels{iChan} = data.channels{1,Input_labels(iChan)};
%         Imp_temp(iChan) = data.Imp_channels(1,Input_labels(iChan));
    end
    data.channels =  temp_channels;
    data.sorted = 'yes';
%     data.Imp_labels = Imp_temp;
end




%%  Create a probe conversion script for the A8x8 probe

if strcmp(probe,'A8x8')
    %     if Connection == 1;
    Input_labels = [37 36 38 35 39 34 40 33 45 44 46 43 47 42 48 41 53 52 54 51 55 50 56 49 61 60 62 59 63 58 64 57 5 4 6 3 7 2 8 1 13 14 16 11 15 10 16 9 21 20 22 19 23 18 24 17 29 28 30 27 31 26 32 25];
    Loc_labels = [5 4 6 3 7 2 8 1 13 14 16 11 15 10 16 9 21 20 22 19 23 18 24 17 29 28 30 27 31 26 32 25 37 36 38 35 39 34 40 33 45 44 46 43 47 42 48 41 53 52 54 51 55 50 56 49 61 60 62 59 63 58 64 57];
    %     elseif Connection ==0;
    %         input_labels = [5 4 6 3 7 2 8 1 13 14 16 11 15 10 16 9 21 20 22 19 23 18 24 17 29 28 30 27 31 26 32 25 37 36 38 35 39 34 40 33 45 44 46 43 47 42 48 41 53 52 54 51 55 50 56 49 61 60 62 59 63 58 64 57];
    %         Loc_labels = [5 4 6 3 7 2 8 1 13 14 16 11 15 10 16 9 21 20 22 19 23 18 24 17 29 28 30 27 31 26 32 25 37 36 38 35 39 34 40 33 45 44 46 43 47 42 48 41 53 52 54 51 55 50 56 49 61 60 62 59 63 58 64 57];
    %     end
    data.labels = Loc_labels;
    temp_channels = cell(1,64);
    for iChan = 1:size(data.channels,2)
        temp_channels{iChan} = data.channels{1,Input_labels(iChan)};
%         Imp_temp(iChan) = data.Imp_channels(1,Input_labels(iChan));
    end
    data.channels =  temp_channels;
    data.sorted = 'yes';
%     data.Imp_labels = Imp_temp;
end

end

