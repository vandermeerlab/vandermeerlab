function [ data_ft, data_filtered ] = AMPX_filter(data, low, high, order)
data_filtered = data;

butter_check % checks to make sure the FT butter function is not being used.  


disp(['Starting Filter between ' num2str(low) '-' num2str(high) ' ...'])
[z,p,k] = butter(order, [low high] * 2 / data.hdr.Fs); % note, we ask for 3 outputs instead of 2
[sos,g] = zp2sos(z,p,k); % convert to SOS format
h = dfilt.df2sos(sos,g); % create filter object
temp_bp_data = cell(1,64);
for iChan = 64:-1:1
    temp_bp_data{1,iChan} = filtfilt(sos,g,data.channels{1,iChan});
end
% put the Bandpassed signal back into the data_ft structure
data_filtered.channels = temp_bp_data; % It is helpful to separate the filtered data from the unfiltered to diagnos any issues in the data. 
clear temp_bp_data;
%% Filter the data into the FieldTrip format
addpath('D:\Users\mvdmlab\My_Documents\GitHub\BIOL680\FieldTrip\external\signal') %This set of functions is needed later on.
data_ft = AMPX_makeft(data_filtered);

disp('Filter Finished');

end

