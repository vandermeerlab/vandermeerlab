function out = ReadNanoZ(cfg_in,fn)
% function out = ReadNanoZ(cfg,fn)
%
% read in NanoZ impedance report
%
% INPUTS:
%  fn: string, name of file to load
%
% OUTPUTS:
%  out.mag: impedance magnitude (in MOhm)
%  out.phase: impedance phase (in degrees)
%
% CONFIGS:
%  cfg.sitemap: integer array, mapping from NanoZ "Site" to in index of out fields 
%   e.g. if sitemap(1) = 18, then NanoZ site 1 values are assigned to
%   out.mag(18) and out.phase(18)
%  cfg.nChannels: integer, number of channels to read
%  cfg.nSkipLines: number of lines to skip at file start before values can
%   be read in
%
% MvdM 2016-02-06

cfg_def = [];
cfg_def.sitemap = [18 27 28 29 17 30 31 32  1  2  3 16  4  5  6 15 ...
                   20 21 22 23 19 24 25 26  7  8  9 14 10 11 12 13];
cfg_def.nSkipLines = 2;
cfg_def.nChannels = 32;
cfg_def.verbose = 1;               

cfg = ProcessConfig(cfg_def,cfg_in);

fid = fopen(fn);

for iL = 1:cfg.nSkipLines
    line = fgetl(fid);
end

for iCh = 1:cfg.nChannels
    line = fgetl(fid);
    
    if line == -1
       error('End of file reached before all channels read.'); 
    end
    
    this_line = sscanf(line,'%d %f %f');
    
    if cfg.verbose
        disp(line);
    end
    
    out.mag(cfg.sitemap(this_line(1))) = this_line(2);
    out.phase(cfg.sitemap(this_line(1))) = this_line(3);

end

fclose(fid);