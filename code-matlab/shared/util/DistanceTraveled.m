function d = DistanceTraveled(cfg_in,pos,iv_in)
% function d = DistanceTraveled(cfg_in,pos,iv_in)
%
% returns distance traveled in pos during interval iv_in
% (sum of absolute diffs), in units of pos
%
% MvdM 2016-10-25


cfg_def = [];
cfg = ProcessConfig2(cfg_def,cfg_in);

nDim = length(pos.label);

% first, construct tsd with dist diffs
switch nDim
    case 1
        diffs = diff(pos.data);
        diffs = tsd(pos.tvec(2:end),abs(diffs));
    case 2
        diffs_x = diff(getd(pos,'x'));
        diffs_y = diff(getd(pos,'y'));
        diffs = tsd(pos.tvec(2:end),sqrt(diffs_x.^2 + diffs_y.^2));
end

these_diffs = restrict(diffs,iv_in);
d = nansum(these_diffs.data);