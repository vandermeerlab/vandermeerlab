function evts_out = AMPX_get_random_evts(data, evts, ExpKeys)
%% AMPX_get_random_evts: samples random events of the same length as the 
%actual gamma events and ensures they are not overlapping. 
%
%   Inputs: 
%     - data: in the AMPX structure
%     - evts: [struct] output from AMPX_Julien_DetectEvents
%
%   Outputs:
%     - evts_out: [struct] contains the original evts with rand_high, and
%     rand_low
%
%% cfg
cfg.f_label = {'low', 'high'};
s = rng; % set the random seed for reproduceability
rng(s)
%%
for iFreq = 1:length(cfg.f_label)
    
r = sort(randi([1 length(data.tvec)], 1, length(evts.(cfg.f_label{iFreq}).tstart)));
for ii = length(r):-1:1
    rand_iv.(cfg.f_label{iFreq}).tstart(ii,1) = data.tvec(r(ii));
    rand_iv.(cfg.f_label{iFreq}).tend(ii,1) = rand_iv.(cfg.f_label{iFreq}).tstart(ii) + evts.(cfg.f_label{iFreq}).tend(ii)-evts.(cfg.f_label{iFreq}).tstart(ii);
end
rand_iv.(cfg.f_label{iFreq}).cfg = evts.(cfg.f_label{iFreq}).cfg;
end



%% compile the random events such that they do not over lap with high or low gamma
ran_iv.tstart = [rand_iv.high.tstart; rand_iv.low.tstart];
ran_iv.tend = [rand_iv.high.tend; rand_iv.low.tend];
ran_iv.cfg = rand_iv.low.cfg;
fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before low Gamma DifferenceIV \n'])

pre_diff = length(ran_iv.tend);
ran_iv = DifferenceIV([], ran_iv, evts.low);
disp([num2str(length(ran_iv.tstart)) '# of ran events after low Gamma DifferenceIV ' ])
post_diff = length(ran_iv.tend);
fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])

fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before high Gamma DifferenceIV \n'])
pre_diff = length(ran_iv.tend);
ran_iv = DifferenceIV([], ran_iv,evts.high);
disp([num2str(length(ran_iv.tstart)) '# of ran events after high Gamma DifferenceIV ' ])
post_diff = length(ran_iv.tend);
fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])

%% loop this until there is only a group of ran_ivs that are the same length as the lg and Hg but not overlapping with either
loop_num = 0;
checks = 10000:1000:400000;
lg_iv = evts.low;
lg_low_iv = evts.low_low_tr;
hg_iv = evts.high;
hg_low_iv = evts.high_low_tr;
cfg.verbose = 0;
while length(ran_iv.tstart) ~= (length(lg_iv.tstart)+length(hg_iv.tstart))
    if ismember(loop_num, checks)
%         hg_iv.tstart = hg_iv.tstart(1:end-5);
%         hg_iv.tend = hg_iv.tend(1:end-5);
%         lg_iv.tstart = lg_iv.tstart(1:end-5);
%         lg_iv.tend = lg_iv.tend(1:end-5);
        s = rng('shuffle');
        disp(['Ran not hitting at ' num2str(loop_num)])
    end
    loop_num = 1+loop_num;
    %     fprintf(['\n\nStarting loop #' num2str(loop_num) '\n\n'])
    
    r = randi([1 length(data.tvec)], 1, ((length(lg_iv.tstart) + length(hg_iv.tstart))- length(ran_iv.tstart)));
    ran_iv_app = [];
    for ii = length(r):-1:1
        ran_iv_app.tstart(ii,1) = data.tvec(r(ii));
        ran_iv_app.tend(ii,1) = ran_iv_app.tstart(ii) + hg_low_iv.tend(ii)-hg_low_iv.tstart(ii);
    end
    ran_iv_app.cfg = ran_iv.cfg;
    ran_iv_app.cfg.seed = s;
    pre_diff = length(ran_iv.tend);
    ran_iv = DifferenceIV(cfg, ran_iv, ran_iv_app);
    post_diff = length(ran_iv.tend);
    %     fprintf(['\n' num2str(pre_diff-post_diff) '  Events removed ran_iv VS ran_iv_app\n'])
    if (pre_diff-post_diff) ==0
        ran_iv.tstart = [ran_iv.tstart; ran_iv_app.tstart];
        ran_iv.tend = [ran_iv.tend; ran_iv_app.tend];
    else
        continue
    end
    %     fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before low Gamma DifferenceIV \n'])
    pre_diff = length(ran_iv.tend);
    ran_iv = DifferenceIV(cfg, ran_iv, lg_low_iv);
    %     disp([num2str(length(ran_iv.tstart)) '# of ran events after low Gamma DifferenceIV ' ])
    post_diff = length(ran_iv.tend);
    %     fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])
    
    %     fprintf(['\n' num2str(length(ran_iv.tstart)) '# of ran events before high Gamma DifferenceIV \n'])
    pre_diff = length(ran_iv.tend);
    ran_iv = DifferenceIV(cfg, ran_iv,hg_low_iv);
    %     disp([num2str(length(ran_iv.tstart)) '# of ran events after high Gamma DifferenceIV ' ])
    post_diff = length(ran_iv.tend);
    %     fprintf([num2str(pre_diff-post_diff) '  Events removed\n'])
    %     disp(['Loop Number ' num2str(loop_num)])
end
fprintf(['\n Finished at loop #' num2str(loop_num) '\n\n'])

%% split the random values into two arrays for low and high gamma
ran_lg_iv = []; ran_hg_iv = [];
r_lg = randperm(length(lg_iv.tstart));
ran_lg_iv.tstart = ran_iv.tstart(r_lg);
ran_lg_iv.tend = ran_iv.tend(r_lg);

ran_hg_iv.tstart = ran_iv.tstart((ismember(ran_iv.tstart, ran_lg_iv.tstart))==0);
ran_hg_iv.tend= ran_iv.tend((ismember(ran_iv.tend, ran_lg_iv.tend))==0);

ran_lg_iv.tstart = sort(ran_lg_iv.tstart);
ran_lg_iv.tend = sort(ran_lg_iv.tend);
ran_hg_iv.tstart = sort(ran_hg_iv.tstart);
ran_hg_iv.tend = sort(ran_hg_iv.tend);


%% collect everything
evts_out = evts;

evts_out.rand_low = ran_lg_iv;
evts_out.rand_high = ran_hg_iv;

evts_out.rand_low.cfg.seed = s;
evts_out.rand_high.cfg.seed = s;

%% report
fprintf(['Low events in: ' num2str(length(evts_out.low.tstart)) '      |    ' num2str(length(evts_out.rand_low.tstart)) ' random "low" epochs'])
fprintf(['High events in: ' num2str(length(evts_out.high.tstart)) '      |    ' num2str(length(evts_out.rand_high.tstart)) ' random "high" epochs'])

end