function ix = findAlignment(D, tstmp, tstmp_unit)

% tsd/findAlignment
% 	ix = findAlignment(D, tstmp)
% 	ix = findAlignment(D, tstmp, tstmp_unit)
%
% 	Returns an index i= such that D(ix) occurs at timestamp tstmp
% 	Finds closest index.
%
%
%INPUTS
%   tstmp:       raw times to find
%   tstmp_unit:  units of tstmp (assumed to be same units as D if not specified)
%
%
% ADR
% version: L5.0
% modified 13 Sep 02 by ncst to use nearest neighbor interpolation
% v5.0: JCJ 3/2/2003 includes support for time units

nIX = length(tstmp);


if isempty(strmatch('units',fieldnames(D)))
    warning('D units not specified: assuming D and tstmp have same units ' )
    ix = interp1(D.t,1:length(D.t),tstmp,'nearest');
    ix = ix(find(~isnan(ix)));
    return
end

if ~exist('tstmp_unit','var')
    tstmp_unit=D.units;
end

if ~strcmp(D.units, tstmp_unit)
    tstmp=convert(tstmp,tstmp_unit,D.units);
end

ix = interp1(D.t,1:length(D.t),tstmp,'nearest');

ix = ix(~isnan(ix));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = convert(R,unit,unitflag)

switch (unitflag)
    case 'sec'
        switch unit
            case 'sec'
                R = R;
            case 'ts'
                R = R/10000;
            case 'ms'
                R = R/1000;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
        
    case 'sec0'
        switch unit
            case 'sec'
                R = (R - min(R));
            case 'ts'
                R = (R - min(R))/10000;
            case 'ms'
                R = (R - min(R))/1000;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
        
    case 'ts'
        switch unit
            case 'sec'
                R = R*10000;
            case 'ts'
                R = R;
            case 'ms'
                R = R*10;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
        
    case 'ms'
        switch unit
            case 'sec'
                R = R*1000;
            case 'ts'
                R = R/10;
            case 'ms'
                R = R;
            otherwise
                warning('tsa has invalid units: no conversion possible' );
        end
    otherwise
        warning('Convert called with invalid unitflag: no conversion possible');
end
