function [values, oldmin, oldmax] = rescale(values, newmin, newmax)
%RESCALE           Rescales a data set.
%   [VALUESOUT, OLDMIN, OLDMAX] = RESCALE(VALUESIN, NEWMIN, NEWMAX), given
%   a VALUESIN array of arbitrary dimension, will return a VALUESOUT of
%   the same size and linearly rescaled to the range NEWMIN:NEWMAX.  The
%   minimum and maximum values of the original matrix are also returned.
%
%   Note: if VALUESIN contains only a single value, the scaled output will
%   be set to (NEWMIN + NEWMAX) / 2.
%
% from the Chronux toolbox

[oldmin oldmax] = minmax(values);

oldrange = oldmax - oldmin;
newrange = newmax - newmin;

% special case
if (oldrange == 0)
    values = repmat((newmax + newmin)/2, size(values));
    return;
end

% Rescaling conceptually uses the following formula:
%     zero_to_one = (oldvalues - oldmin) ./ oldrange;
%     newvalues   = (zero_to_one .* newrange) + newmin;
% But doing it out longhand like this takes two matrix multiplications
% and two matrix additions.  The command below takes a more convoluted
% route to do the same thing with half the number of full matrix operations.
scale = oldrange ./ newrange;
shift = oldmin - (scale .* newmin);
values = (values - shift) ./ scale;


function [xmn,xmx,mni,mxi] = minmax(x)
%CORE_MINMAX       Core computational routine for MINMAX.
%   [XMN,XMX,MNI,MXI] = CORE_MINMAX(X) returns scalars such that
%   XMN = X(MNI) = min(X) and XMX = X(MXI) = max(X).  Ties in indices are
%   broken in favor of the lowest magnitude index.
%   
%   NaN values are ignored unless the input is all NaN.  In this case, XMN
%   and XMX are set to NaN, while MNI and MXI are set to 1.  This mimics
%   the behavior of the Matlab native MIN/MAX functions. 
%
%   CONDITIONS
%   ----------
%   X must be a real vector of type DOUBLE.  An N-D array X is treated as X(:).
%   Infinite values are allowed.   NaN's are ignored (see above).

[xmn,mni] = min(x(:));
[xmx,mxi] = max(x(:));

return;