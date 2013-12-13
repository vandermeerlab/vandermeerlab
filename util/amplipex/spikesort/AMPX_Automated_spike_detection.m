%%%%%% Automated rough spike detection %%%%%%%
% This script will go through all of the recording sessions specified in
% the 'sessions_to_detect.m' set of strings which are in the format
% 'yyyy-mm-dd-[session type]' (eg: 2013-10-22-pre).  It will load each
% session and then bases on the session type it will process different
% components of the session.  For the pre and post sessions it will use the
% entire session.  The main session will use both the last 5 mins and the
% final omisson trial.  It will take pictures of each group of 8 channels
% for the desired epoch in each session condition as well as the
% correlation matrix for the data to help with referencing.


%% set the path and any variables
addpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab')
addpath('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\util\amplipex\spikesort')
cd('D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\util\amplipex\spikesort')

%% Create the session list if not already done

% for ii = 1:length(sessions)
%     if ii <6
%     sessions_to_detect{ii} = ['2013-08-0' num2str(ii+4)];
%     else
%        sessions_to_detect{ii} = ['2013-08-' num2str(ii+4)];
%     end
% end


%% load the sessions to be detected.
rat_id = 36;
fname = strrep('R0name_sessions_to_detect','name',num2str(rat_id));
load(fname,'-mat')
dir_name = ['R0' num2str(rat_id)];



for folder_ind = 3:length(sessions)  % loops thorugh all the folders in the session list
    current_folder = [dir_name '-' sessions{folder_ind}];
    cd(['D:\DATA\R0' num2str(rat_id) '\' current_folder]);
    % detect the avalability of the main data set recording
    if exist([current_folder '.dat']) ~=0
        sess_info.data = 'yes';
    else
        sess_info.data = 'no';
    end
    % detect the avalability of the pre record
    if exist([current_folder '-pre.dat']) ~=0
        sess_info.pre = 'yes';
    else
        sess_info.pre = 'no';
    end
    
    % detect the avalability of the post recorde
    if exist([current_folder '-post.dat']) ~=0
        sess_info.post = 'yes';
    else
        sess_info.post = 'no';
    end
    
    if strcmp(sess_info.pre, 'yes')==1;
        AMPX_spike_images(current_folder,'type','pre')
    end
    if strcmp(sess_info.post, 'yes')==1;
        AMPX_spike_images(current_folder,'type','post')
    end
    if strcmp(sess_info.data, 'yes')==1;
        AMPX_spike_images(current_folder,'type','main')
    end

    
    %%%%%%%%%% load the pre record and apply the detection function (to be made
    %%%%%%%%%% later)

        
        
    
end  % end the loop thourgh all the folders in the session  list.


