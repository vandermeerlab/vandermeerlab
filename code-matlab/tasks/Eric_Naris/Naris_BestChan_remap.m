function [best, all_best] = Naris_BestChan_remap(ExpKeys, varargin)
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
            best_chans = [64 63 56 55 48 62  54 47 61 ];
        elseif strcmp(location, 'dm')
            best_chans = [1 9 2 10 17 3 11 18 3];
        elseif strcmp(location, 'dl')
            best_chans = [57 58 49 50 59 41 51 60 42];
        elseif strcmp(location, 'vm')
            best_chans = [8 16 7 15  24 6 14 23 22];
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
            best_chans = [64 56 63 55 62 48 54 47 61 ];
        elseif strcmp(location, 'dm')
            best_chans = [1 9 2 10 17 3 11 18 3];
        elseif strcmp(location, 'dl')
            best_chans = [57 58 49 50 59 41 51 60 42];
        elseif strcmp(location, 'vm')
            best_chans = [8 16 7 15  24 6 14 23 22];
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
