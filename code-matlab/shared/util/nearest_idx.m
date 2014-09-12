function yi = nearest_idx(x,y)
% function yi = nearest_idx(x,y)
%
% returns indices yi such that x-y(yi) is minimized
%
% NOTE: x and y MUST BE SORTED! NO WARNING GIVEN FOR SPEED!
%
% MvdM 2014-06-25 adapted from:
% http://stackoverflow.com/questions/2142826/mapping-2-vectors-help-to-vectorize
% 

yi = zeros(length(x), 1);
last_M = 1;
for N = 1:length(x)
  % search through M until we find a match.
  for M = last_M:length(y)
    dist_to_curr = abs(y(M) - x(N));
    dist_to_next = abs(y(M+1) - x(N));

    if dist_to_next > dist_to_curr
      yi(N) = M;
      last_M = M;
      break
    else
      continue
    end

  end % M
end % N

%% fix
%if M<numel(xm) 
%    dist_to_next = abs(xm(M+1) - xn(N)); 
%else xmap4(N) = M; 
%    break 
%end
%dist_to_curr = abs(xn(N) - xm(M));