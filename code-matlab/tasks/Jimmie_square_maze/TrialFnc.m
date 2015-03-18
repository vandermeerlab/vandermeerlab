function TrialInfo = TrialFnc(eventLog)
% eventLog to TrialInfo

clear TrialInfo

cue_type = 1; %manual input for now, modified control script to record it in eventLogs
for i = 1:eventLog.trials_initiated
    %     TrialInfo.cue_type(i) = eventLog.cue_type(i);
    switch cue_type
        case 1
            TrialInfo.cue_type = 'light';
            TrialInfo.cue_ID(i,1) = eventLog.light_ID(i);            
        case 2
            TrialInfo.cue_type = 'sound';
            TrialInfo.cue_ID(i,1) = eventLog.sound_ID(i);            
    end
    TrialInfo.rewarded(i,1) = eventLog.rewarded(i);
    if TrialInfo.rewarded(i,1) == 2;
        TrialInfo.rewarded(i,1) = 1;
    end
    % TrialInfo.approach(i,1) = eventLog.approach(i);
    TrialInfo.trialT(i,1) = eventLog.trialT(i);
    
    temp = find (eventLog.pb_breaksT == eventLog.trialT(i),1,'first'); %find the corresponding time in eventLog.pb_breaksT of the trial initiation
    if eventLog.pb_breaksID(temp + 1) == eventLog.pb_breaksID(temp) + 4 %see if the next pb break is the adjacent reward pb
        temp2 = find (eventLog.ALLnosepokeT == eventLog.pb_breaksT(temp +1),1,'first'); %find the location of this nosepoke in the eventLog.ALLnosepokeT variable
        TrialInfo.nosepokeT(i,1) = eventLog.pb_breaksT(temp + 1);
        TrialInfo.trial_to_nosepokeT(i,1) = TrialInfo.nosepokeT(i,1) - TrialInfo.trialT(i,1);
        for jj = temp2:length(eventLog.unnosepokeT)
            if eventLog.ALLnosepokeID(temp2) ~= eventLog.unnosepokeID(jj) %find the last unnosepoke from current reward pb
                TrialInfo.unnosepokeT(i,1) = eventLog.unnosepokeT(jj - 1);
                TrialInfo.nosepoke_length(i,1) = TrialInfo.unnosepokeT(i,1) - TrialInfo.nosepokeT(i,1);
                if i ~= eventLog.trials_initiated
                    TrialInfo.unnosepoke_to_trialT(i,1) = eventLog.trialT(i + 1) - TrialInfo.unnosepokeT(i,1);
                else
                    TrialInfo.unnosepoke_to_trialT(i,1) = 0;
                end
                temp3 = eventLog.unnosepokeT(jj - 1) - eventLog.ALLnosepokeT(temp2); %find time of the last unnosepoke
                break 
            end
        end
        if temp3 > 1
            TrialInfo.approached(i,1) = 1;
            TrialInfo.quick_approach(i,1) = 0;
        else
            TrialInfo.approached(i,1) = 0;
            TrialInfo.quick_approach(i,1) = 1;
        end
    else
        TrialInfo.nosepokeT(i,1) = 0;
        TrialInfo.unnosepokeT(i,1) = 0;
        TrialInfo.trial_to_nosepokeT(i,1) = 0;
        TrialInfo.nosepoke_length(i,1) = 0;
        TrialInfo.approached(i,1) = 0; 
        TrialInfo.quick_approach(i,1) = 0;
        TrialInfo.unnosepoke_to_trialT(i,1) = 0;
    end
    
    if i ~= eventLog.trials_initiated
        TrialInfo.trial_length(i,1) = eventLog.trialT(i + 1) - eventLog.trialT(i);
        temp4 = find (eventLog.backwardsT > TrialInfo.trialT(i,1),1,'first');
        if eventLog.backwardsT(temp4) < eventLog.trialT(i + 1);
            TrialInfo.backwards(i,1) = 1;
        else
            TrialInfo.backwards(i,1) = 0;
        end
    else
        TrialInfo.trial_length(i,1) = eventLog.ALLnosepokeT(end) - eventLog.trialT(i);
        TrialInfo.backwards(i,1) = 0;
    end
    
    TrialInfo.summary(i,1) = TrialInfo.cue_ID(i,1);
    TrialInfo.summary(i,2) = TrialInfo.rewarded(i,1);
    TrialInfo.summary(i,3) = TrialInfo.approached(i,1);
    TrialInfo.summary(i,4) = TrialInfo.trialT(i,1);
    TrialInfo.summary(i,5) = TrialInfo.nosepokeT(i,1);
    TrialInfo.summary(i,6) = TrialInfo.unnosepokeT(i,1);
    TrialInfo.summary(i,7) = TrialInfo.trial_length(i,1);
    TrialInfo.summary(i,8) = TrialInfo.trial_to_nosepokeT(i,1);
    TrialInfo.summary(i,9) = TrialInfo.nosepoke_length(i,1);
    TrialInfo.summary(i,10) = TrialInfo.unnosepoke_to_trialT(i,1);
    TrialInfo.summary(i,11) = TrialInfo.quick_approach(i,1);
    TrialInfo.summary(i,12) = TrialInfo.backwards(i,1);
    TrialInfo.summary_legend{1,1} = 'cue ID';
    TrialInfo.summary_legend{1,2} = 'rewarded';
    TrialInfo.summary_legend{1,3} = 'approached';
    TrialInfo.summary_legend{1,4} = 'trial start time';
    TrialInfo.summary_legend{1,5} = 'nosepoke time';
    TrialInfo.summary_legend{1,6} = 'unnosepoke time';
    TrialInfo.summary_legend{1,7} = 'trial to trial length';
    TrialInfo.summary_legend{1,8} = 'trial to nosepoke length';
    TrialInfo.summary_legend{1,9} = 'nosepoke length';
    TrialInfo.summary_legend{1,10} = 'unnosepoke to trial length';
    TrialInfo.summary_legend{1,11} = 'quick approach';
    TrialInfo.summary_legend{1,12} = 'backwards running';
end