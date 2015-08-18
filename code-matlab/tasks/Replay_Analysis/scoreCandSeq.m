function [score] = scoreCandSeq(cfg_in,CAND_iv_in,S)
%% SCORECANDSEQ Score Candidate Seuences
% SCORE = scoreCandSeq(cfg_in,CAND_iv_in,S,pos) returns a score of the candidate replay
% event sequences. SCORE is a struct containing different possible scoring methods (i.e.,
% p-values).
%
%   INPUT:
%       cfg_in - input cfg parameters
%       CAND_iv_in - input candidate event iv
%       S - input spiketrain
%       pos - input position vector
%
%   OUTPUT:              
%       score - struct containing replay scores
%
%   CFG OPTIONS:
%       cfg.method = 'ROCorr'; %'RODist'
%       cfg.nShuffles = 300; %number of shuffles for bootstrapping
%       cfg.seqMethod = 'mean'; %('exact','first','mean') determines which spikes to use for sequencing
%       cfg.cdfcheck = 0;
%       cfg.display = 0;
%       cfg.siglvl = 0.05; % alpha for significance testing
%
% youkitan 2014-12-20

%% Parse input cfg parameters

cfg_def.method = 'ROCorr'; %'RODist'
cfg_def.nShuffles = 250; %number of shuffles for bootstrapping
cfg_def.seqMethod = 'mean'; %('exact','first','mean') determines which spikes to use for sequencing
cfg_def.cdfcheck = 0;
cfg_def.display = 0;
cfg_def.siglvl = 0.05; % alpha for significance testing
cfg = ProcessConfig2(cfg_def,cfg_in);
 
%% Get sequence orderings and templates for spikes during CAND_iv
sequences = {};
for iEvent = 1:length(CAND_iv_in.tstart)
    S_event = restrict(S,CAND_iv_in.tstart(iEvent),CAND_iv_in.tend(iEvent));
    spk = [];
    inds = [];

    for iC = 1:length(S_event.t)
        if ~isempty(S_event.t{iC})
            switch cfg.seqMethod
                case 'exact' %use all observed spikes
                    spk = cat(1,spk,S_event.t{iC}); 
                    inds = cat(1,inds,iC*ones(1,numel(S_event.t{iC}))');
                    
                case 'first' %time of first spike
                    spk = cat(1,spk,S_event.t{iC}(1));
                    inds = cat(1,inds,iC);        
        
                case 'mean' %average spike time (peak of tc?)
                    spk = cat(1,spk,mean(S_event.t{iC}));
                    inds = cat(1,inds,iC);        
            end %sequence type
        end %empties
    end %iterate cells
    
    [~,idx] = sort(spk);

    sequences{1,iEvent} = inds(idx); % spike sequence (observed ordering)
    sequences{2,iEvent} = inds; % template sequence (place field ordering)
end %iterate events

% Create shuffled sequences and calculate rank order correlations
rhos_shuffled = zeros(cfg.nShuffles,length(sequences));
pvals_shuffled = rhos_shuffled;
rhos = zeros(1,length(sequences));
pvals = rhos;

for iSeq = 1:length(sequences)
    s1 = sequences{1,iSeq};
    s2 = sequences{2,iSeq};
    if isempty(s1) || isempty(s2)
        s1 = nan; s2 = nan;
    end
    [rhos(iSeq),pvals(iSeq)] = corr(s1,s2,'type','Spearman');
    
    if iSeq==1 || mod(iSeq,10)==0
        str = sprintf('Calculating correlations for sequence %d',iSeq);
        disp(str);
    end

    
    for iSh = 1:cfg.nShuffles
        rand_idx = randperm(length(sequences{2,iSeq}));
        s2 = sequences{2,iSeq}(rand_idx);
        [rhos_shuffled(iSh,iSeq),pvals_shuffled(iSh,iSeq)] = corr(s1,s2,'type','Spearman');
    end
end

if strcmp(cfg.method,'ROCorr')
    keep_idx = zeros(1,length(rhos));
    
    for iEvt = 1:length(rhos)
        
        xvals = -1:0.05:1;
        rhos_shuffled_dist = rhos_shuffled(:,iEvt);
        [n1,x1] = hist(rhos_shuffled_dist,xvals);

        % Calculate alpha values for two tails
        a_neg = prctile(rhos_shuffled_dist,2.5); %<2.5%
        a_pos = prctile(rhos_shuffled_dist,97.5); %>97.5%       
        
        if rhos(iEvt) < a_neg || a_pos < rhos(iEvt)
            significance = 'yes';
        else
            significance = 'no';
        end
        
        if pvals(iEvt) < cfg.siglvl && strcmp(significance,'no')
            sprintf('conflicting significance')
        elseif pvals(iEvt) > cfg.siglvl && strcmp(significance,'yes')
            sprintf('conflicting significance')
        end
            
            
            
        if cfg.display
            % Compare shuffled distribution to observed value
            bar(x1,n1,'FaceColor',[.7 .7 .7]); hold on;
            ylims = get(gca,'Ylim');
            plot([rhos(iEvt) rhos(iEvt)],[ylims(1) ylims(2)],'r');
            plot([a_neg a_neg],[ylims(1) ylims(2)],':g');
            plot([a_pos a_pos],[ylims(1) ylims(2)],':g');
            axis([-1.1 1.1 ylims]);

            % Annotate each event
            str1 = sprintf('Rho-value: %4f',rhos(iEvt));
            str2 = sprintf('p-value: %4f',pvals(iEvt));
            str3 = sprintf('Signigicant (vs shuffle): %s',significance);
            text(-1,ylims(2)-5,{str1,str2,str3},'FontSize',8);
            title(sprintf('Event %d', iEvt));
            waitforbuttonpress
            cla;
        end
        
        % Add significant events to keep_idx
        if strcmp(significance,'yes')
            keep_idx(iEvt) = 1;	
        end
    end
    close

    score.ROcorr = rhos;
    score.pvals = pvals;
    score.significance = logical(keep_idx);
    
else
    % Plot distributions
    xvals = -1:0.1:1;
    rhos_shuffled_dist = rhos_shuffled(:);
    [n1,x1] = hist(rhos_shuffled_dist,xvals);
    [n2,x2] = hist(rhos,xvals);

    figure;
    bar(x1,n1,'FaceColor',[.7 .7 .7]); hold on;
    bar(x2,n2,'r');
    axis([-1.1 1.1 get(gca,'Ylim')]);

    % Two-sample Kolmogorov-Smirnov test
    [hypothesis,pvalue] = kstest2(rhos,rhos_shuffled_dist);

    % Store replay distribution score
    score.hypothesis = hypothesis;
    score.pvalue = pvalue;

end
%% check CDFs to make sure K-S test is working
if cfg.cdfcheck
    figure;
    cdfplot(rhos); hold all;
    cdfplot(rhos_shuffled_dist);

    x = -1:0.001:1;
    sigma = 0.44;
    norm = normcdf(x,0,sigma);
    plot(x,norm);

    legendnames = {'data','shuffled',['normal dist (\sigma = ' num2str(sigma) ')']};
    legend(legendnames,'Location','Best')
    hold off;
end        



end