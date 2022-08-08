%% visualize full model predictions vs. data
cfg.LFPrange = [-2 1];
cfg.spikeRange = [0 0.1];
cfg.predictionRange = 1; % scaling factor
cfg.predictionErrorRange = 1;

iC = 1;

% LFP
csc_resc = resc(csc.data,cfg.LFPrange);
plot(csc.tvec,csc_resc,'Color',[0.7 0.7 0.7]);
hold on;

% spike train
plot(TVECc(te_idx),resc(spk_binned(te_idx),cfg.spikeRange),'.g','MarkerSize',10);

% prediction
plot(TVECc(te_idx),m3p,'.k','MarkerSize',10);

% prediction error
plot(TVECc(te_idx),sse3,'.r','MarkerSize',10);

%% model difference

%%
function out = resc(in,range)

if length(range) == 1 % scaling factor
    out = in .* range;
elseif length(range) == 2 % range
    % first scale to [0 1]
    out = in - min(in);
    out = out ./ max(out);
    
    % then use input
    rdiff = diff(range);
    out = out .* rdiff;
    out = out + range(1);
    
else
   error('Unknown input length.'); 
end




end