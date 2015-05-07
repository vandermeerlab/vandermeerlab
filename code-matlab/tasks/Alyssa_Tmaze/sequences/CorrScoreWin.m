function [score] = CorrScoreWin(cfg_in,iv_in,S)
% function [score] = CorrScoreWin(cfg_in,iv_in,S)
% returns a score of the candidate replay event sequences. SCORE is a struct containing different possible scoring methods (i.e.,
% p-values).
%
%   INPUT:
%       cfg_in - input cfg parameters
%       iv_in - input candidate event iv
%       S - input spiketrain
%
%   OUTPUT:              
%       score - struct containing replay scores
%
%   CFG OPTIONS:'
%       cfg.nShuffles = 300; %number of shuffles for bootstrapping
%       cfg.seqMethod = 'mean'; %('exact','first','mean') determines which spikes to use for sequencing
%       cfg_def.shuffleType = 3; % as per ShuffleTS()
%       cfg_def.dt = 0.025;
%       cfg_def.twin = 0.025:0.025:0.15;
%
% youkitan 2014-12-20, MvdM 2015

%% Parse input cfg parameters

cfg_def.nShuffles = 250; %number of shuffles for bootstrapping
cfg_def.seqMethod = 'exact'; %('exact','first','mean') determines which spikes to use for sequencing
cfg_def.shuffleType = 3; % as per ShuffleTS()
cfg_def.dt = 0.025;
cfg_def.twin = 0.025:0.025:0.2;
cfg_def.nMinCells = 5;

cfg = ProcessConfig2(cfg_def,cfg_in);
 
%% Get sequence orderings and templates for spikes during CAND_iv
nEvt = length(iv_in.tstart);

WIN_rho_perc = nan(nEvt,1); 
WIN_rho_obs = nan(nEvt,1); 
WIN_rho_obs_pval = nan(nEvt,1);
WIN_start_idx = nan(nEvt,1);
WIN_end_idx = nan(nEvt,1);
WIN_iv = iv([0 0]);
WIN_rho_shuf = nan(nEvt,cfg.nShuffles);
SHUF_full_rho_perc = {};
OBS_rho = {};
OBS_pval = {};
OBS_tstart = {};
OBS_tend = {};

parfor iEvt = 1:nEvt
    
    fprintf('Event %d/%d...\n',iEvt,nEvt);
    
    % original S
    full_S = restrict(S,iv_in.tstart(iEvt),iv_in.tend(iEvt));
    
    % create shuffled S's so that we don't have to do this for each
    % start-end combination later
    cfg_shuf = []; 
    cfg_shuf.t0 = iv_in.tstart(iEvt); cfg_shuf.t1 = iv_in.tend(iEvt);
    cfg_shuf.mode = cfg.shuffleType;
    
    shuf_S = cell(cfg.nShuffles,1);
    for iSh = 1:cfg.nShuffles
        shuf_S{iSh} = ShuffleTS(cfg_shuf,full_S);
    end
    
    
    starts = iv_in.tstart(iEvt):cfg.dt:iv_in.tend(iEvt);
    
    ends = cat(2,starts(1) + cfg.twin,iv_in.tend(iEvt));
    ends = ends(ends <= iv_in.tend(iEvt));
    
    % init vars for this start-end pair
    SE_rho = nan(length(starts),length(ends)); SE_rho_perc = nan(length(starts),length(ends));
    SE_pval = nan(length(starts),length(ends)); 
    
    SE_rho_shuf = nan(cfg.nShuffles,length(starts),length(ends));
    SE_pval_shuf = nan(cfg.nShuffles,length(starts),length(ends));
    
    for iStart = 1:length(starts)
        
        this_start = starts(iStart);

        for iEnd = 1:length(ends)
            
            this_end = ends(iEnd);
            if this_end <= this_start, continue; end
            
            this_S = restrict(full_S,this_start,this_end);
            
            [OBS_seq_temp,OBS_templ_temp] = GetSeq(cfg,this_S);
            
            % NOTE could use this conditional to exclude insufficient
            % neurons -- unique() on template
            if length(unique(OBS_seq_temp)) < cfg.nMinCells % not enough spikes in interval!
               SE_rho_perc(iStart,iEnd) = 0; % worst possible percentile
               SE_rho(iStart,iEnd) = 0;
               SE_pval(iStart,iEnd) = Inf;
               continue; 
            end
            
            [SE_rho(iStart,iEnd),SE_pval(iStart,iEnd)] = corr(OBS_seq_temp,OBS_templ_temp,'type','Spearman');           
            
            if isnan(SE_rho(iStart,iEnd)) % happens if all ranks are the same, for instance
                SE_rho_perc(iStart,iEnd) = 0; % worst possible percentile
                SE_rho(iStart,iEnd) = 0;
                SE_pval(iStart,iEnd) = Inf;
                continue;
            end
            
            SHUF_rho_temp = nan(cfg.nShuffles,1); SHUF_pval_temp = nan(cfg.nShuffles,1);
            
            for iSh = 1:cfg.nShuffles
                
                temp_S = restrict(shuf_S{iSh},this_start,this_end);
                                                                
                [SHUF_seq,SHUF_templ] = GetSeq(cfg,temp_S);
                
                if length(SHUF_seq) <= 2 % not enough spikes in interval!
                    SHUF_rho_temp(iSh) = 0; SHUF_pval_temp(iSh) = Inf;
                    continue;
                end

                %[SHUF_rho_temp(iSh),SHUF_pval_temp(iSh)] = corr(SHUF_seq,SHUF_templ,'type','Spearman');
                SHUF_seq = tiedrank(SHUF_seq); SHUF_templ = tiedrank(SHUF_templ);
                SHUF_rho_temp(iSh) = fastcorr(SHUF_seq,SHUF_templ);
                
            end
            
            % get percentile of current start-end pair
            SE_rho_perc(iStart,iEnd) = nansum(abs(SHUF_rho_temp) < abs(SE_rho(iStart,iEnd)))./cfg.nShuffles;

            % store all thee shuffled values and p-values (to keep track of
            % scored shuffles)
            SE_rho_shuf(:,iStart,iEnd) = SHUF_rho_temp; SE_pval_shuf(:,iStart,iEnd) = SHUF_pval_temp;
            
        end % of ends
    end % of starts
    
    % now store what needs to be kept for this evt
    SHUF_full_rho_perc{iEvt} = SE_rho_perc;
    
    OBS_rho{iEvt} = SE_rho;
    OBS_pval{iEvt} = SE_pval;
    OBS_tstart{iEvt} = starts;
    OBS_tend{iEvt} = ends;
    
    % note, could use some better way to handle ties
    [rho_val,rho_idx] = max(SE_rho_perc(:)); [rho_x,rho_y] = ind2sub(size(SE_rho_perc),rho_idx);
    
    WIN_rho_perc(iEvt) = rho_val;
    WIN_rho_shuf(iEvt,:) = sq(SE_rho_shuf(:,rho_x,rho_y));
    WIN_pval_shuf(iEvt,:) = sq(SE_pval_shuf(:,rho_x,rho_y));
    WIN_rho_obs(iEvt) = SE_rho(rho_x,rho_y);
    WIN_rho_obs_pval(iEvt) = SE_pval(rho_x,rho_y);
    WIN_iv(iEvt) = iv([starts(rho_x) ends(rho_y)]);
    WIN_start_idx(iEvt) = rho_x; WIN_end_idx(iEvt) = rho_y;
    
end %iterate events

score.SHUF_full_rho_perc = SHUF_full_rho_perc;
score.OBS_rho = OBS_rho;
score.OBS_pval = OBS_pval;
score.OBS_tstart = OBS_tstart;
score.OBS_tend = OBS_tend;
score.WIN_rho_obs = WIN_rho_obs;
score.WIN_rho_perc = WIN_rho_perc;
score.WIN_rho_obs_pval = WIN_rho_obs_pval;
score.WIN_iv = WIN_iv;
score.WIN_start_idx = WIN_start_idx;
score.WIN_end_idx = WIN_end_idx;
score.WIN_rho_shuf = WIN_rho_shuf;
score.WIN_pval_shuf = WIN_pval_shuf;

function [seq,templ] = GetSeq(cfg,S_event)

spk = []; % could speed up by knowing how many spikes there are going to be
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
                
            case 'peak'
                error('Peak method not yet implemented');
                
        end % sequence type
    end % empties
end % iterate cells

[~,idx] = sort(spk);

seq = inds(idx); % spike sequence (observed ordering)
templ = inds; % template sequence (place field ordering)

function Rho = fastcorr(x,y)

% http://www.mathworks.com/matlabcentral/newsreader/view_thread/119195
m = size(x,1);
xc = x - repmat(sum(x,1)/m,m,1); % remove mean
yc = y - repmat(sum(y,1)/m,m,1);
Sigma = (xc' * yc) / (m-1);
sx = sqrt(sum(xc.^2,1) / (m-1)); % compute std dev
sy = sqrt(sum(yc.^2,1) / (m-1));
Rho = Sigma ./ (sx'*sy);