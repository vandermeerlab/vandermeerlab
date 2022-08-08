function out = CalcWVDistances(cfg_in, out, out2)
% out = CalcWVDistances(cfg_in, out, out2)
%
% compute distances between all cell waveforms in out and all cell waveforms in out2
% (out and out2 are the output of CategorizeStriatumWave.m)

cfg_def = [];
cfg = ProcessConfig(cfg_def, cfg_in);

n1 = size(out.wv, 1);
n2 = size(out2.wv, 1);

for iX = 1:n1
    for iY = 1:n2
        
        wv1 = out.wv(iX, :);
        wv2 = out2.wv(iY, :);
        
        out.SSE(iX, iY) = sum((wv1 - wv2).^2); % sum of squares
        out.SSEn(iX, iY) = out.SSE(iX, iY) ./ sum((wv1 + wv2).^2); % normalized sum of squares
        out.abs_diff(iX, iY) = sum(abs(wv1 - wv2)); % summed absolute difference
        out.abs_diffn(iX, iY) = out.abs_diff(iX, iY) ./ sum((wv1 + wv2)); % normalized absolute difference
                
        temp = corrcoef(wv1, wv2);
        out.corr(iX, iY) = temp(1, 2); % correlation
        
        [peak1_val, peak1_idx] = max(wv1);
        out.peakdiff(iX, iY) = peak1_val - wv2(peak1_idx); % raw difference compared to peak of wv1
        out.peakdiffn(iX, iY) = out.peakdiff(iX, iY) ./ peak1_val; % normalized difference (wv2 is what fraction of wv1 peak)
        
    end
end
