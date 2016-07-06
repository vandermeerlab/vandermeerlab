function AMPX_Naris_generate_figures()
%% AMPX_Naris_generate_figures: Generates all the figures used in Carmichael et al. 
%   
%     - figure 1
%         - C: raw traces of a random low gamma event
%         - D: heat plot of the gamma power during the same random low gamma event
%         - E: 4 x 2 of the average high and low gamma power across
%              sessions per rat (1-4)
%         - S1: same as E but the color axis scale is the min/max across
%              all sessions/rats
%
%     - figure 2 [AMPX_Naris_fig_2_task.m]
%         - A: generates figure 2 which contains the heatmaps
%              for the average rewarded/approach gamma 50 and gamma 80 
%              events for the two rats that reached task criterion
%         - B: plane fitting histograms and statistics (R^2 and
%              angles/directions)
%
%     - figure 3 [AMPX_Naris_fig_3_phase.m]
%         - A: creates 4 x 2 figues of the phase differences
%              for both across entire gamma events 
%         - B: 4 x 2 for phase differences across the central three cycles.
%     
%     - figure 4 [AMPX_Naris_fig_4_CSD.m]
%         - A: kernal current source density for the example event used in
%              fig 1E
%         - B: 4 x 2 of the average kCSD for the central three cycles in
%              each event per session per rat
%
%     - figure 5 [AMPX_Naris_fig_5_gamma_count.m]
%         - bar plot of the normalized count of gamma events across all
%           Naris occlusion recordings
%         - figure
%
%
%EC - 2016-06-14

%% load the pre and post data:

load('C:\temp\Naris_all_data_pre.mat');

load('C:\temp\Naris_all_data_post.mat');


%% Figure 1: individual gamma example (raw and heat map), averages across all rats

AMPX_Naris_fig_1_example(all_data_pre, all_data_post, 'Example', 29, 'session_name', 'R061_2014_09_26');


%% Figure 2: 2x2x2 task

AMPX_Naris_fig_2_task(all_data_pre, all_data_post)%, 'save_fig', 'yes');


%% Figure 3: creates 4 x 2 figues of the phase differences for both across entire gamma events and within the triplets

AMPX_Naris_fig_3_phase([], all_data_pre, all_data_post, 'save_fig', 'yes');

%% Figure 4: CSD for same event used in figure 1 C/D and average for each rat (50/80)

AMPX_Naris_fig_4_CSD([], 'Example', 29, 'session_name', 'R061_2014_09_26')

%% Figure 5: Naris experiments: 

% three example sessions "white" PSDs. 


% Power comparisons


% Gamma count bar graph
AMPX_Naris_fig_5_gamma_count()
