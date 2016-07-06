%% Naris Paper Sandbox

%% list of sessions to analyze
% Session_list = {'R054-2014-10-11', 'R054-2014-10-12', 'R054-2014-10-13', 'R054-2014-10-14'};
Session_list = {'R054-2014-10-13', 'R054-2014-10-14', 'R054-2014-10-15', ...
    'R049-2014-02-07', 'R049-2014-02-08', 'R049-2014-02-10',... % 'R049-2014-02-09',...   
    'R061-2014-09-26', 'R061-2014-09-27', 'R061-2014-09-28',...
    'R045-2014-04-16', 'R045-2014-04-17', 'R045-2014-04-18'};

type = 'pre';

%% list of sessions to analyze
% Session_list = {'R054-2014-10-11', 'R054-2014-10-12', 'R054-2014-10-13', 'R054-2014-10-14'};
Session_list = {'R049-2014-02-07', 'R049-2014-02-08', 'R049-2014-02-10',... % 'R049-2014-02-09',...   
    'R045-2014-04-16', 'R045-2014-04-17', 'R045-2014-04-18'};

type = 'post';

%% get the task sessions
Session_list = {'R049-2014-02-07', 'R049-2014-02-08', 'R049-2014-02-10',... 
    'R045-2014-04-16', 'R045-2014-04-17', 'R045-2014-04-18'};
type = 'task';
%% loop over each session to get: events, power, event_phase, middle cycles, cycle_phase

for iSess =1 %:length(Session_list)
    disp(['Running Session: ' (strrep(Session_list{iSess},'-', '_'))])
    if strcmp(type, 'pre') == 1
        all_data.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'post') == 1
        all_data_post.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    elseif strcmp(type, 'task') ==1
       all_data_task.(strrep(Session_list{iSess},'-', '_')) = AMPX_Naris_pipeline([Session_list{iSess}], 'session_type', type, 'plane_plot', 'yes');
    end
end
 
%% save the all_data struct
    if strcmp(type, 'pre') == 1
        save('C:\temp\Naris_all_data_pre.mat', 'all_data_pre', '-v7.3')
    elseif strcmp(type, 'task') ==1
        save('C:\temp\Naris_all_data_task.mat', 'all_data_task', '-v7.3')
    elseif strcmp(type, 'post') ==1
        save('C:\temp\Naris_all_data_post.mat', 'all_data_post', '-v7.3')
    end
%% get R^2 across sessions
sess_list = fieldnames(all_data_pre);
bands = {'lg','lg_ran', 'hg', 'hg_ran'};% 'hg'};
all.lg.rsq = []; all.hg.rsq = [];all.lg_ran.rsq = []; all.hg_ran.rsq = [];
all.lg.rsq2 = []; all.hg.rsq2 = [];all.lg_ran.rsq2 = []; all.hg_ran.rsq2 = [];
for iSess = 1:length(sess_list)
    for iband = 1:length(bands)
        all.(bands{iband}).rsq = [all.(bands{iband}).rsq, all_data_pre.(sess_list{iSess}).(bands{iband}).power.stats.rsq];
        all.(bands{iband}).rsq2 = [all.(bands{iband}).rsq2, all_data_pre.(sess_list{iSess}).(bands{iband}).power.stats.rsq2];
    end
    hist_lg = AMPX_plane_count_hist(all_data_pre.(sess_list{iSess}).lg.power.stats, all_data_pre.(sess_list{iSess}).lg_ran.power.stats, 'low');
    hist_lg = AMPX_plane_count_hist(all_data_pre.(sess_list{iSess}).hg.power.stats, all_data_pre.(sess_list{iSess}).hg_ran.power.stats, 'high');
    disp(sess_list{iSess})
pause; close all

end

hist_lg = AMPX_plane_count_hist(all.lg, all.lg_ran, 'low');
hist_hg = AMPX_plane_count_hist(all.hg, all.hg_ran, 'high');

%% power figures

h_power_evt = AMPX_Naris_fig_1_example([], all_data, 'R061-2014-09-27', 'Example', 110)

%% phase figure

AMPX_Naris_fig_3_phase([], all_data)
