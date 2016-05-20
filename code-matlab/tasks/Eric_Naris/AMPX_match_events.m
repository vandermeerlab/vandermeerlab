function ctrl_evts = AMPX_match_events(cfg_in, evts, detect_csc)
%% AMPX_match_events: wrapper function for MatchGammaEvents which:
%     " finds event times with low gamma power matching gamma events" see
%     MatchGammaEvents help
%
%
%  INPUTS:
%     - cfg [struct] contains parameters for control event width
%     - evts [struct]: output from AMPX_JulienDetectEvents
%
%  Outputs
%     - ctrl_evts [struct]: ivs for the outputs from MatchGammaEvents

%% set up default parameters
cfg_def.PARAM_twin = 0.2; % half-width of time window
cfg_def.PARAM_control_dt = 5; % offset for control classifier
cfg_def.debug = 0;
cfg = ProcessConfig(cfg_def, cfg_in);

%% run the low gamma power matching epochs


% create fixed-size events
lg_t = IVcenters(evts.low); hg_t = IVcenters(evts.high);
this_lg = iv(lg_t-cfg.PARAM_twin,lg_t+cfg.PARAM_twin); %this_lg_control = iv(lg_t-PARAM_twin+PARAM_control_dt,lg_t+PARAM_twin+PARAM_control_dt);
this_hg = iv(hg_t-cfg.PARAM_twin,hg_t+cfg.PARAM_twin); %this_hg_control = iv(hg_t-PARAM_twin+PARAM_control_dt,hg_t+PARAM_twin+PARAM_control_dt);

% find matched events, as a control
cfg_match = []; cfg_match.twin = [-cfg.PARAM_control_dt cfg.PARAM_control_dt]; cfg_match.evt_twin = [-cfg.PARAM_twin cfg.PARAM_twin];
ctrl_evt = MatchGammaEvents(cfg_match,detect_csc,lg_t,hg_t);
this_lg_control = iv(ctrl_evt.lg-cfg.PARAM_twin,ctrl_evt.lg+cfg.PARAM_twin);
this_hg_control = iv(ctrl_evt.hg-cfg.PARAM_twin,ctrl_evt.hg+cfg.PARAM_twin);


ctrl_evts.lg = this_lg_control;
ctrl_evts.hg = this_hg_control;

%%
if cfg.debug == 1
    cfg_plot = []; cfg_plot.display = 'tsd';
    PlotTSDfromIV(cfg_plot,this_lg,detect_csc);
    cfg_plot.iv_only = 1; cfg_plot.fgcol = 'g';
    PlotTSDfromIV(cfg_plot,this_lg_control,detect_csc);
end


%% print the output
fprintf(['\nlg events in:           ' num2str(length(evts.low.tstart)) '\n'])
fprintf(['Control low events out: ' num2str(length(ctrl_evts.lg.tstart)) '\n'])
fprintf(['\nhg events in:            ' num2str(length(evts.high.tstart)) '\n'])
fprintf(['Control high events out: ' num2str(length(ctrl_evts.hg.tstart)) '\n'])