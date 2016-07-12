function laps = get_LOOCV_laps(cfg_in,available_laps,this_lap,this_EncLap)
% function laps = get_LOOCV_laps(cfg_in,available_laps,this_lap,this_EncLap)
%
% return this_EncLap (this is a count) idxs taken from available_laps (a
% set of idxs), excluding this_lap (a single idx)
%
% example:
%
% current lap is 8, need 5 closest laps from set of 10
%
% >> get_LOOCV_laps([],1:10,8,5)
%
% ans =
%
%     9     7    10     6     5
%
% MvdM 2016-07-01

cfg_def = [];
cfg_def.mode = 'nearest'; % 'nearest' or 'from_end' or 'dist'

cfg = ProcessConfig(cfg_def,cfg_in);

available_laps = sort(available_laps);
this_lap_idx = find(available_laps == this_lap);

if isempty(this_lap_idx) 
    error('current lap not included in available laps!');
end

switch cfg.mode
    
    case 'from_end'
        if this_EncLap > length(available_laps)
            error('Not enough laps available!');
        end
        
        % just count from the back
        available_laps = setxor(available_laps,this_lap);
        laps = available_laps(end:-1:end-(this_EncLap-1));
        
    case 'nearest'
        if this_EncLap > length(available_laps)
            error('Not enough laps available!');
        end
        
        % alternative: add laps closest in distance
        available_laps = setxor(available_laps,this_lap);
        
        lap_dist = abs(available_laps-this_lap);
        lap_dist(available_laps > this_lap) = lap_dist(available_laps > this_lap) - 0.1; % symmetry breaking to ensure later laps are cnosen first when equal distance
        
        [~,sorted_lap_idx] = sort(lap_dist);
        
        keep_idx = sorted_lap_idx(1:this_EncLap);
        laps = available_laps(keep_idx);

    case 'dist'
        % return lap specified distance apart, if available
        if this_EncLap > length(available_laps)-1
            disp('Not enough laps available!');
            laps = []; return;
        end

        lap_dist = available_laps-this_lap;
        keep_idx = find(lap_dist == this_EncLap);
        laps = available_laps(keep_idx);
        
end