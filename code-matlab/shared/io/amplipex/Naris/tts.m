function tetrodes = tts(cfg);
% functional purpose: to call up the channels for requested tetrode
tetrode = cfg.tts_to_process;
rID = cfg.fname(1:4);
if strcmp(rID, 'R053')
list_of_channels{1} = [22 24 26 28];
list_of_channels{2} = [21 23 25 27];
list_of_channels{3} = [14 16 18 20];
list_of_channels{4} = [13 15 17 19];
list_of_channels{5} = [6 8 10 12];
list_of_channels{6} = [5 7 9 11];
list_of_channels{7} = [1 2 3 4];
list_of_channels{8} = [67 68 69 70]; %fill space
list_of_channels{9} = [54 56 58 60];
list_of_channels{10} = [53 55 57 59];
list_of_channels{11} = [46 48 50 52];
list_of_channels{12} = [45 47 49 51];
list_of_channels{13} = [38 40 42 44];
list_of_channels{14} = [37 39 41 43];
list_of_channels{15} = [30 32 34 36];
list_of_channels{16} = [29 31 33 35];

tetrodes.channels_array = [];

for iT = 1:1:length(tetrode)
    tetrodes.identity{iT} = tetrode(iT);
    tetrodes.channels{iT} = list_of_channels{tetrode(iT)};
    tetrodes.channels_array = cat(2, tetrodes.channels_array, tetrodes.channels{iT});
end
elseif strcmp(rID, 'R060')
    list_of_channels{1} = [26 28 30 32];
list_of_channels{2} = [25 27 29 31];
list_of_channels{3} = [18 20 22 24];
list_of_channels{4} = [17 19 21 23];
list_of_channels{5} = [10 12 14 16];
list_of_channels{6} = [9 11 13 15];
list_of_channels{7} = [2 4 6 8];
list_of_channels{8} = [1 3 5 7];
list_of_channels{9} = [58 60 62 64];
list_of_channels{10} = [57 59 61 63];
list_of_channels{11} = [50 52 54 56];
list_of_channels{12} = [49 51 53 55];
list_of_channels{13} = [42 44 46 48];
list_of_channels{14} = [41 43 45 47];
list_of_channels{15} = [34 36 38 40];
list_of_channels{16} = [33 35 37 39];

tetrodes.channels_array = [];

for iT = 1:1:length(tetrode)
    tetrodes.identity{iT} = tetrode(iT);
    tetrodes.channels{iT} = list_of_channels{tetrode(iT)};
    tetrodes.channels_array = cat(2, tetrodes.channels_array, tetrodes.channels{iT});
end
end
