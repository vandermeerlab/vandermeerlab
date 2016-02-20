function eql_iv = GetEqualBehaviorIV(cfg_in,metadata)
% function eql_iv = GetEqualBehaviorIV(cfg_in,metadata)
%
% returns iv with only those rest intervals preceded by "equal behavior":
% same number of L and R trials
%
% CONFIGS:
% cfg_def = [];
% cfg_def.excludeBlockedTrials = 0;
%
% MvdM 2015-11-16 initial version

%
cfg_def = [];
cfg_def.excludeBlockedTrials = 0;

cfg = ProcessConfig(cfg_def,cfg_in);

%
nTrials = length(metadata.taskvars.sequence);

tstart = []; tend = [];
Neql_iv = 0; % counter for number of rest_ivs with equal behavior

for iT = 1:nTrials
    
    if EqualBehavior(metadata.taskvars.sequence,iT)
        
        % before adding, check that this rest iv indeed comes after preceding trial
        if ~metadata.taskvars.rest_iv.tstart(iT) > metadata.taskvars.trial_iv.tend(iT)
            error('Rest iv start appears before trial end');
        end
        
        Neql_iv = Neql_iv + 1;
        tstart(Neql_iv) = metadata.taskvars.rest_iv.tstart(iT);
        tend(Neql_iv) = metadata.taskvars.rest_iv.tend(iT);
        
    end % of equal behavior match
    
end

eql_iv = iv(tstart,tend);

function iseql = EqualBehavior(sequence,nTrial)

this_sequence = sequence(1:nTrial);

nL = sum(strcmp(this_sequence,'L'));
nR = sum(strcmp(this_sequence,'R'));

iseql = nL == nR;