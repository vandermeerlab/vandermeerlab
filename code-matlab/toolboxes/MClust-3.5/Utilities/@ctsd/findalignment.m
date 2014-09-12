function ix = findAlignment(D, tstmp, tstmp_unit)

% tsArray/findAlignment
% 	ix = findAlignment(D, tstmp)
% 	ix = findAlignment(D, tstmp, tstmp_unit)
%
% 	Returns an index i= such that D(ix) occurs 
%		at timestamp tstmp
%
%INPUTS
%   tstmp:      raw times to find
%   tstmp_unit: units of tstmp (assumed to be same units as D if not specified)
%
%
% ADR
% version: L5.0
% v5.0: JCJ 2/27/2003 includes support for time units


if isempty(strmatch('units',fieldnames(D)))
    warning('D units not specified: assuming D and tstmp have same units ' )
    ix = floor((tstmp - D.t0)/D.dt)+1;
    ix(find(ix==0))=1;  % if ix ever equals 0 then there will be errors further on
    return
end

if ~exist('tstmp_unit','var')
    tstmp_unit=D.units;
end

if strcmp(D.units, tstmp_unit)
    ix = floor((tstmp - D.t0)/D.dt)+1;
else
    T0=starttime(D,tstmp_unit);
    DT=dt(D,tstmp_unit);
    
    ix = floor((tstmp - T0  )/ DT )+1;
end

ix(ix>length(D.data))=length(D.data);  % if ix ever larger than the highest index in D then there will be errors further on
ix(ix<1)             =1;               % if ix ever equals 0 then there will be errors further on


