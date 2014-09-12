function spk_t = lastSpike(S)
% function spk_t = lastSpike(S)
%
% returns time of last spike

spk = vertcat(S.t{:});
spk_t = max(spk);