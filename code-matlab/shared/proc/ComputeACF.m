function [acf,tvec] = ComputeACF(cfg,spk_binned)

if isfield(cfg,'maxlag')
    [acf,tvec] = xcorr(spk_binned,spk_binned,cfg.maxlag);
else
    [acf,tvec] = xcorr(spk_binned,spk_binned);
end

tvec = tvec.*cfg.binsize;
acf(ceil(length(acf)/2)) = 0;

if isfield(cfg,'sided')
    cut_idx = ceil(length(acf)/2);
    switch cfg.sided
        case 'one' % only uses data from center to end.  
            acf = acf(cut_idx:end);
            tvec = tvec(cut_idx:end);
        case 'onezero' % keeps the length for convolution
            acf(1:cut_idx-1) = 0;
    end
end

end