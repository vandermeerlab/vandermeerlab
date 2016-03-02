function nearestvalues = nearval3(queryvalues, lookupvalues, direction)
%nearval3   Returns the values in lookupvalues that are closest to 
%  the numbers in queryvalues. 
%
%  nearestvalues = nearval3(queryvalues, lookupvalues, direction)
%
%  INPUTS
%  queryvalues: the value(s) to look for in lookupvalues
%  lookupvalues: list of ordered values (ascending)
%  direction: 
%   0: default, closest values
%   1: look forwards in lookupvalues
%  -1: look backwards in lookupvalues   
%  If a value falls outside the range of lookupvalues, it always chooses 
%  the closest value, regardless of direction.
%
%  ACarey, Dec 2014
%  Relies on nearest_idx3 by MvdM


%% 

switch nargin
    
    case 2
        nearestvalues = lookupvalues(nearest_idx3(queryvalues,lookupvalues,0));
        
    case 3
        nearestvalues = lookupvalues(nearest_idx3(queryvalues,lookupvalues,direction));
end
        
end

