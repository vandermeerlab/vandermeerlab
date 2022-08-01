function tvec = read_smi(fin_name)
% load timestamps in .smi files saved by Cheetah when recording video mp4
% files
%
% MvdM 2021

cfg_def = [];
cfg_def.initNframes = 1e6;

fid = fopen(fin_name,'rt');

fout_name = [];
line = 1;
frame_no = 1;

tvec = nan(cfg_def.initNframes, 1);

while line ~= -1
    
    line = fgetl(fid);
    
    if line == -1
        break;
    end
    
    match = regexp(line, 'SCC>\d*<', 'match');
    if ~isempty(match)
       tvec(frame_no) = str2double(match{1}(5:end-1)); 
       frame_no = frame_no + 1;
    end
        
end

tvec = tvec(1:frame_no-1);
tvec = tvec*1e-10; % this may be wrong

fprintf('%d frames loaded.\n', frame_no);