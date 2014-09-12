function spk_t = firstSpike(S)
% function spk_t = firstSpike(S)
%
% returns time of first spike

spk = vertcat(S.t{:});
spk_t = min(spk);