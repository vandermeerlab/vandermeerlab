function [ ExpKeys ] = emi_loadExpKeys(unique_folder)
% Loads expKeys from data folder. Returns ExpKeys object.
% Usage:
% unique_folder = 'R068-2014-12-01_recording';
% ExpKeys = emi_loadExpKeys(unique_folder);

key_id = strcat(unique_folder(1:4),'_',unique_folder(6:9),'_',...
    unique_folder(11:12),'_',unique_folder(14:15));

unique_key = strcat(key_id,'_ExpKeys.m');
run(unique_key);

save([unique_folder(1:15),'-ExpKeys.mat'],'ExpKeys');
end