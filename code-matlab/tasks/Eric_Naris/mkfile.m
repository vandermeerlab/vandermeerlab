function [file_name] = mkfile()
%% make_file_name

current_dir = cd;
file_name = current_dir(end-14:end);
clear current_dir
end
