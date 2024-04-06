function [tsd_out, bin] = resampleTSD(cfg_in, tsd_in, new_tvec)
% function [tsd_out, bin] = resampleTSD(cfg_in, tsd_in, new_tvec)
%
% resamples tsd_in using time points defined in new_tvec
%
% outputs:
%  tsd_out: resampled output TSD
%  bin: [1 x tsd_in.tvec] vector containing output bin assignment of
%  each input sample (rebinning mode only)
% 
% inputs:
%  tsd_in: original input TSD
%  new_tvec: [1 x nTimePoints] column vector with new sample times
%
% options:
%
% 'rebinning' method averages data points in tsd_in that are nearest 
%   to the times in new_tvec to make data points in tsd_out
%
% 'interp' method simply uses MATLAB interp function
%
% cfg_def.method = 'rebinning'; % 'rebinning', 'interp', 'resample'
% cfg_def.interp_method = 'nearest'; % uses interp options
%
% MvdM 2024

cfg_def = [];
cfg_def.method = 'rebinning'; % 'rebinning', 'interp', 'resample'
cfg_def.interp_method = 'nearest'; % uses interp options

cfg = ProcessConfig(cfg_def, cfg_in);

if ~CheckTSD(tsd_in)
    error('Input not a correctly formed tsd');
end

if ~iscolumn(new_tvec)
    warning('new_tvec input is not a column vector, converting...')
    new_tvec = new_tvec';
end

nCh = size(tsd_in.data, 1);

tsd_out = tsd_in; tsd_out.tvec = new_tvec; tsd_out.data = [];

switch cfg.method

    case 'rebinning'

        % first, set up bin edges
        tvec_diff = diff(new_tvec);
        tvec_edges = [-Inf; new_tvec(1:end-1) + tvec_diff/2; Inf];

        % now find mapping from input samples to output bins
        [counts, ~, bin] = histcounts(tsd_in.tvec, tvec_edges);

        nBins = max(bin);
        if nBins > length(new_tvec), error('DEBUG: Too many output bins (%d)', nBins); end

        for iCh = nCh:-1:1

            for iBin = nBins:-1:1 % slow, how to speed up?
                tsd_out.data(iCh,iBin) = nanmean(tsd_in.data(iCh, bin == iBin));
            end % bins

        end % channels

        % check if any new bins have no samples
        if any(counts == 0)
            warning('%d resampled time points have no data', sum(counts == 0));
        end

        % check if any new bins have NaNs
        if any(isnan(counts))
            warning('NaNs in resampled data');
        end

    case 'interp'

        for iCh = nCh:-1:1

            tsd_out.data(iCh,:) = interp1(tsd_in.tvec, tsd_in.data(iCh,:), new_tvec, cfg.interp_method);

        end % channels

    case 'resample'

        error('Not yet implemented')

end % methods


tsd_out = History(tsd_out, mfilename, cfg);
