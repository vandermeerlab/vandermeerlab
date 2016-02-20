function values = rescmean(values,newmean)
%RESCMEAN Rescale data according to a new mean
%   VALUESOUT = rescmean(VALUESIN,NEWMEAN) given a VALUESIN array of 
%   arbitrary dimension, will return a VALUESOUT of the same size and 
%   linearly rescaled such than the mean sits at NEWMEAN. 
%
%   VALUESIN can also be a tsd; in this case, VALUESIN.data is rescaled.
%
% A.Carey, Feb 2015.

if isfield(values,'data')
    if mean(values.data) == 0
        warning('Initial mean is zero, cannot rescale')
    else
        values.data = values.data .* newmean/mean(values.data);
    end
    
else
    if mean(values) == 0
        warning('Initial mean is zero, cannot rescale')
    else
        values = values .* newmean/mean(values);
    end
    
end

