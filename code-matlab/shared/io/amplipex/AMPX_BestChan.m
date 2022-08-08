function [best, all_best] = AMPX_BestChan(ExpKeys, varargin)
if strcmp(ExpKeys.ProbeType, 'A8x8') ==0;
    error('This is only working for A8x8 at this time')
end

%% define the varibales
location = 'all';
all_loc{1} = 'dl'; all_loc{2} = 'dm'; all_loc{3} = 'vl'; all_loc{4} = 'vm';

extract_varargin

%%
if strcmp(location, 'all')
    for iL = 1:length(all_loc)
        location = all_loc{iL};
        
        if strcmp(location, 'vl')
            best_chans = [57 64 49 56 58 50 15 63 4 9 55 59 2 6];
        elseif strcmp(location, 'dm')
            best_chans = [38 35 46 37 43 27 36 45 25 18 40 44 29 20 14];
        elseif strcmp(location, 'dl')
            best_chans = [61 53 60 52 7 62 59 54 5 14 51 11 63 16];
        elseif strcmp(location, 'vm')
            best_chans = [34 33 39 47 17 42 41 30 28 40 19 48 36 28];
        end
        
        %%
        best_chan = [];
        for iChan  = 1:length(best_chans)
            if isempty(intersect(best_chans(iChan), ExpKeys.BadChannels)) ==1
                best_chan = [best_chan;  best_chans(iChan)];
            else
                continue
            end
        end
        best = best_chan(1);
        all_best.(all_loc{iL}) = best;
    end
else
    if strcmp(location, 'vl')
        best_chans = [57 64 49 56 58 50 15 63 4 9 55 59 2 6];
    elseif strcmp(location, 'dm')
        best_chans = [38 35 46 37 43 27 36 45 25 18 40 44 29 20 14];
    elseif strcmp(location, 'dl')
        best_chans = [61 53 60 52 7 62 59 54 5 14 51 11 63 16];
    elseif strcmp(location, 'vm')
        best_chans = [34 33 39 47 17 42 41 30 28 40 19 48 36 28];
    end
    
    %%
    best_chan = [];
    for iChan  = 1:length(best_chans)
        if isempty(intersect(best_chans(iChan), ExpKeys.BadChannels)) ==1
            best_chan = [best_chan;  best_chans(iChan)];
        else
            continue
        end
    end
    best = best_chan(1);
    all_best.(location) = best;    
end
