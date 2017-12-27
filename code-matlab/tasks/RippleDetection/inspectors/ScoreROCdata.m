function [bestF1Score,bestJScore] = ScoreROCdata(cfg,ROCdata)
%SCOREROCDATA Compute F1-score and Youden J-score of the ROCs.
%
%   [bestF1Score,bestJScore] = ScoreROCdata(cfg,ROCdata)
%
%   CONFIG OPTIONS
%      cfg.use1s = 1;
%      cfg.use2s = 1;
%      cfg.use3s = 1;
%      cfg.use4s = 0;
%      cfg.use5s = 0;
%          1 - Consider annotated intervals of a given rating (1,2,3,4,5)
%          0 - Do not consider.
%
%     cfg.verbose = 1; 
%          1 - Print informative text to the command window 
%          0 - Don't print text to the command window
%
%   OUTPUTS
%   
%   bestF1Score is the highest possible F1 Score over all thresholds
%   tested. We consider only 1/2/3/4 SWRs and ignore the 5s. F1 score is a
%   measure of test accuracy that values false positives and negatives
%   equally. F-score of 1 is a perfect test, F-score of 0 is a random
%   guesser, F-score of -1 is a pathological failure.
%
%   See https://en.wikipedia.org/wiki/F1_score
%
%   bestJScore is the highest possible Youden J statistic over all
%   thresholds tested. Again, only considering the 1/2/3/4s. This is just
%   the sum of sensitivity and specificity - 1. i.e. it is the maximum
%   vertical distance from the line y=x on the ROC charts, where the
%   1/2/3/4s have all been properly combined into a single category.
%   Again, 1 is optimal and 0 is a completely random guesser.
%
%   See https://en.wikipedia.org/wiki/Youden%27s_J_statistic
%
% Elyot Grant Dec 2017 initial version
% aacarey Dec 2017

% Parse cfg parameters
    cfg_def.verbose = 1;
    cfg_def.use1s = 1;
    cfg_def.use2s = 1;
    cfg_def.use3s = 1;
    cfg_def.use4s = 0;
    cfg_def.use5s = 0;
    
    mfun = mfilename;
    cfg = ProcessConfig(cfg_def,cfg,mfun);

    % Check cfg.use inputs
    assert(cfg.use1s == 1 || cfg.use1s == 0)
    assert(cfg.use2s == 1 || cfg.use2s == 0)
    assert(cfg.use3s == 1 || cfg.use3s == 0)
    assert(cfg.use4s == 1 || cfg.use4s == 0)
    assert(cfg.use5s == 1 || cfg.use5s == 0)
    
    % Do the thing
        
    numThresholds = length(ROCdata.hitrates1);
    
    bestJScore = -1;
    bestF1Score = -1;
    bestF2Score = -1;
    bestF05Score = -1;
    bestJThr = NaN;
    bestF1Thr = NaN;
    bestF2Thr = NaN;
    bestF05Thr = NaN;
    
    for iThr = 1:numThresholds
        totalTruePositives = cfg.use1s * ROCdata.nEvt1(iThr)+...
                             cfg.use2s * ROCdata.nEvt2(iThr)+...
                             cfg.use3s * ROCdata.nEvt3(iThr)+...
                             cfg.use4s * ROCdata.nEvt4(iThr)+...
                             cfg.use5s * ROCdata.nEvt5(iThr);
        totalFalsePositives = ROCdata.nFalsePos(iThr);
        totalReportedTrue = cfg.use1s * ROCdata.nEvt1(iThr)/ROCdata.hitrates1(iThr) +...
                            cfg.use2s * ROCdata.nEvt2(iThr)/ROCdata.hitrates2(iThr) +...
                            cfg.use3s * ROCdata.nEvt3(iThr)/ROCdata.hitrates3(iThr) +...
                            cfg.use4s * ROCdata.nEvt4(iThr)/ROCdata.hitrates4(iThr) +...
                            cfg.use5s * ROCdata.nEvt5(iThr)/ROCdata.hitrates5(iThr);
        totalFalseNegatives = totalReportedTrue - totalTruePositives;
        
        precision = totalTruePositives/(totalTruePositives + totalFalsePositives);
        recall = totalTruePositives/(totalTruePositives + totalFalseNegatives);
        
        beta = 1; %F1 score
        f1Score = (1+beta^2)*(precision*recall)/(precision*beta^2 + recall);
        beta = 2; %F1 score
        f2Score = (1+beta^2)*(precision*recall)/(precision*beta^2 + recall);
        beta = 0.5; %F1 score
        f05Score = (1+beta^2)*(precision*recall)/(precision*beta^2 + recall);
        
        sensitivity = recall;
        specificity = 1 - ROCdata.falseposrates(iThr);
        jScore = sensitivity + specificity - 1;
        
        if (f1Score > bestF1Score)
            bestF1Score = f1Score;
            bestF1Thr = iThr;
        end
        
        if (f2Score > bestF2Score)
            bestF2Score = f2Score;
            bestF2Thr = iThr;
        end
        
        if (f05Score > bestF05Score)
            bestF05Score = f05Score;
            bestF05Thr = iThr;
        end
        
        if (jScore > bestJScore)
            bestJScore = jScore;
            bestJThr = iThr;
        end
    end
    
    % Get thresholds from config history
    cfg_temp = []; cfg_temp.verbose = 0; cfg_temp.target = 'GetROCdata'; cfg_temp.parameter = 'thresholds';
    thresholds = GetHistory(cfg_temp,ROCdata);
    
    bestF1Thr = thresholds{1,1}(bestF1Thr);
    bestF2Thr = thresholds{1,1}(bestF2Thr);
    bestF05Thr = thresholds{1,1}(bestF05Thr);
    bestJThr = thresholds{1,1}(bestJThr);
    
    % talk to me or not
    if cfg.verbose
        fprintf('%s: best F1 score is %0.4f with a threshold of %0.2f; best J score is %0.4f with a threshold of %0.2f\n',mfun,bestF1Score,bestF1Thr,bestJScore,bestJThr)
        fprintf('%s: best F2 score is %0.4f with a threshold of %0.2f; best F0.5 score is %0.4f with a threshold of %0.2f\n',mfun,bestF2Score,bestF2Thr,bestF05Score,bestF05Thr)
    end
end

