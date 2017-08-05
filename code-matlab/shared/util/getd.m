function [out,idx] = getd(in_tsd,target_field)
% function [out,idx] = getd(in_tsd,target_field)
%
% tsd/ts helper function to return data corresponding to label target_field
%
% MvdM 2014-06-24

idx = find(strcmp(target_field,in_tsd.label));

if isempty(idx)
    error('target field not found.');
else
    if isfield(in_tsd,'data')
        out = in_tsd.data(idx,:);
    elseif isfield(in_tsd,'t') % ts
        out = in_tsd.t{idx};
    else
       error('unknown data format.'); 
    end
end