function [cfg] = Naris_cfgs(cfg)
% cfg.fname = mkfile;
if strcmp(cfg.fname(1:4), 'R053')
    cfg.tetrodes_to_process = 4;
    cfg.chan = 1;
elseif strcmp(cfg.fname(1:4),'R060')
    cfg.tetrodes_to_process = 15;
    cfg.piri_tetrodes_to_process = 6;
    cfg.chan = 4;
elseif  strcmp(cfg.fname(1:4),'R065')
    cfg.tetrodes_to_process = 4;
    cfg.piri_tetrodes_to_process = 6;
    cfg.chan = 1;
end
