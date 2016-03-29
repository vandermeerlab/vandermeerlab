%% make_file_name_pre

current_dir = cd;
file_name = current_dir(end-14:end);
file_name = [file_name '-post'];
clear current_dir
