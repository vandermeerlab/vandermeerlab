function [TSD,oldmin,oldmax] = NormalizeTSD(cfg,TSD,type,newval)
%NORMALIZETSD Rescale or normalize TSD.data
%
%   TSD = NormalizeTSD(cfg,TSD,'mean',newmean)
%   TSD = NormalizeTSD(cfg,TSD,'median',newmedian)
%   TSD = NormalizeTSD(cfg,TSD,'range',[newmin newmax])
%         Rescale TSD.data to fall between newmin and newmax
%
%   For 'mean' and 'median' we take abs(TSD.data) before normalizing, so it
%   doesn't make sense to specify a newmean or newmedian that is negative.
%
%   CFG OPTIONS     
%       cfg.verbose = 1 (default)
%          1 - Print informative text to the command window 
%          0 - Don't print text to the command window
%
% aacarey Feb 2015 initial version of rescmean
% aacarey Dec 2017 renamed RescaleTSD to NormalizeTSD and combined with rescmean

% The inputs (type and newval) don't make sense as cfg defaults...because
% what is a good default? Instead, they are required inputs to the function
cfg_def.verbose = 1;
mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg,mfun);


switch type
    case 'range'
        assert(isnumeric(newval) && length(newval) == 2)
        [TSD.data,oldmin,oldmax] = rescale(TSD.data,newval(1),newval(2));
    case 'median'
        assert(isnumeric(newval) && length(newval) == 1)
        oldmin = min(TSD.data); oldmax = max(TSD.data);
        
        if median(TSD.data) == 0
            % cannot divide by zero, do nothing
            warning('Initial median is zero, cannot normalize')
        else
            TSD.data = TSD.data .* newval/median(abs(TSD.data));
        end
    case 'mean'
        assert(isnumeric(newval) && length(newval) == 1)
        oldmin = min(TSD.data); oldmax = max(TSD.data);
        
        if mean(TSD.data) == 0
            % cannot divide by zero, do nothing
            warning('Initial mean is zero, cannot normalize')
        else
            TSD.data = TSD.data .* newval/mean(abs(TSD.data));
        end
        
    otherwise
        error('Unrecognized type specified')
end

% the inputs aren't configs, but save them in history anyway
cfg.type = type;
cfg.newval = newval;
TSD = History(TSD,mfun,cfg);

end

