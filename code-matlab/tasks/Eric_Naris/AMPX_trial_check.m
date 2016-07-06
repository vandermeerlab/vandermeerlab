function AMPX_trial_check(data_ft)
%% plots 9 random traisl to check the trialification
%
%    INPUTS:
%      - data_ft [struct] trialified data in the Ft format
%


%% make an plot with 9 random trials

ids = randperm(length(data_ft.trial));
if length(ids) >= 16    
ids = ids(1:16);
end
if findall(0, 'type', 'figure') == 200; close(200); end
figure(200)
maximize
for ii  = 1:length(ids)
    subplot(4,4,ii)
    plot(data_ft.time{ids(ii)}, data_ft.trial{ids(ii)}(64,:), 'linewidth', 2, 'color', 'k')
    axis off
end
