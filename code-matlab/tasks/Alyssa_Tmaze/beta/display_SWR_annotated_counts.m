
% Update the data folder and the filename of the annotation ('manualIV')
% file, then run the code to display the count of each category of
% manually-annotated SWR

%% To update
data_folder = 'C:\shared\DATA\Example_mvdmlab_data\R064\R064-2015-04-18';
this_ann_fn = 'R064-2015-04-18-manualIV_post_nofilter_ED_annoted.mat';

%% Main code
cur_dir = pwd();
% go to the data directory
cd(data_folder);
% load the data
evt = load(this_ann_fn);
evt = evt.evt;
% convert to numbers and count the events in each category
ann_mat = str2num(cell2mat(evt.usr.annotation));
categories = sort(unique(ann_mat));
[categories_counts, edges] = histcounts(ann_mat);
% convert in table for ease of plotting & display its contents
SWR_count_table = table(categories, categories_counts', ...
    VariableNames= {'SWR categories', 'Count'})
fprintf('Total events found: %d \n', sum(SWR_count_table.Count));

% return to the original folder
cd(cur_dir);
