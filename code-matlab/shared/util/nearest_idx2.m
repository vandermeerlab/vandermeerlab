function yi = nearest_idx2(x,y)
% function yi = nearest_idx2(x,y)
%
% returns indices yi such that x-y(yi) is minimized
%
% does not require sorted x and y

yi = nan(size(x));

for ix = 1:length(x)
    
   d = abs(y-x(ix));
   [~,yi(ix)] = min(d);
    
end