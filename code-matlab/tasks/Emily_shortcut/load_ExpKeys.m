function [ ExpKeys ] = load_ExpKeys(rat_id,unique_folder)
% Loads expKeys from data folder. Returns ExpKeys object.
% Usage:
% rat_id = 'R068_EI';
% unique_folder = 'R068-2014-12-01_recording';
% ExpKeys = emi_loadExpKeys(rat_id,unique_folder);

cd(fullfile('H:\data-working\Shortcut',rat_id,unique_folder));
key_id = strcat(unique_folder(1:4),'_',unique_folder(14:15));

unique_key = strcat('emi_',key_id,'_expKeys.m');
run(unique_key);

end