function [ outlier_index ] = position_constraint( position_data,window_size,outlier )

outlier_index = [];

for left=1:(length(position_data)-window_size)
    mean_window = mean(position_data(left:left+window_size));
    if abs(position_data(left)-mean_window) > outlier
        outlier_index(end+1) = left;
    end
end
end