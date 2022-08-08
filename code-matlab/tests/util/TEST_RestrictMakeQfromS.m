%% test commutativity of Restrict and MakeQfromS
one_ts_data = (0:1:100)'; % very regular spike times
S = ts; S.t{1} = one_ts_data; S.label{1} = 'label1';

% set up trials
trials = iv;
trials.tstart = 1;
trials.tend = 2;
%trials.tstart = [1 5];
%trials.tend = [2 6];

% Restrict all trials first, then MakeQfromS
Sr = restrict(S, trials);
Qr = MakeQfromS([], Sr);

% Restrict each trial first, then MakeQfromS
Qr2 = tsd;
for iT = 1:length(trials.tstart)
    
   this_S = restrict(S, trials.tstart(iT), trials.tend(iT));
   this_Q = MakeQfromS([], this_S);   
   
   Qr2 = UnionTSD([], Qr2, this_Q);
   
end

% MakeQfromS first, then Restrict
Qr3 = MakeQfromS([], S);
Qr3 = restrict(Qr3, trials);

assert(isequal(length(Qr.tvec), length(Qr2.tvec), length(Qr3.tvec)));