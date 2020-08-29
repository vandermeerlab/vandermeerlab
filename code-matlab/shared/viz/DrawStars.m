function DrawStars(cfg_in, p, coord)
% function DrawStars(cfg, p, coord)

cfg_def.fs = 12;

cfg = ProcessConfig(cfg_def, cfg_in);

if p < 0.001
    txt = '***';
elseif p < 0.01
    txt = '**';
elseif p < 0.05
    txt = '*';
else
    txt = '';
end

th = text(coord(1), coord(2), txt, 'FontSize', cfg.fs);
set(th, 'HorizontalAlignment', 'center')
      