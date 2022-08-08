function [predicted_x, lambda_x, x_binned] = MakeTC_1D(cfg, tvec, data, TVECc, spk_binned)
% construct tuning curve (mean value of spk_binned for binned values of
% data
%
% INPUTS:
%
% tvec, data: describe tuning variable (e.g. position, phase)
% TVECc, spk_binned: describe spike train
%
% cfg.interp: option to use for interpolation, 'linear', 'nearest' (for
%   things like phase)
% cfg.bins: bin edges (vector) for constructing tuning curves
%
% OUTPUTS:
%
% x_binned: binned version of data input
% lambda_x: tuning curve (size cfg.bins), gives spike probability as
%   function of binned tuning variable
% predicted_x: predicted value of tuning variable for each TVECc input
%   sample (simply looked up in lambda_x)

% get variable of interest on common timebase 
x_interp = interp1(tvec, data, TVECc, cfg.interp);

% bin variable of interest
[x_occ_hist, ~, x_binned] = histcounts(x_interp, cfg.bins); 
x_binned(x_interp < cfg.bins(1) | x_interp >= cfg.bins(end)) = NaN;

% edges for binned version of tuning variable; minus 1 to get centers 
% from previous histcounts use, + 0.5 for adding edge
bin_edges = 0.5:1:(length(cfg.bins) - 1) + 0.5; 

if islogical(spk_binned) % binned spike train (binary), can use fast method
    x_spk = x_binned(spk_binned); % variable of interest at spike times
    x_spk_hist = histcounts(x_spk, bin_edges);
    lambda_x = x_spk_hist ./ x_occ_hist;
else % continuous input, need to average across samples in each bin (slow)
    for iB = (length(cfg.bins) - 1):-1:1 % minus 1 to get centers, not edges
        this_samples = (x_binned == iB);
        lambda_x(iB) = nanmean(spk_binned(this_samples));
    end
end

% based on tuning curve, get lambdas for each time bin
predicted_x = nan(size(x_binned));
keep = x_binned ~= 0 & ~isnan(x_binned);
predicted_x(keep) = lambda_x(x_binned(keep));

end