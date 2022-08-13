function [T, WV] = LoadTT_UMS2K(fn, records_to_get, records_flag)
% MClust 4+ loader for UMS2K output (spk struct) saved as .mat file
% MvdM 22

if nargin == 1
    load(fn);
    T = double(spk.unwrapped_times);
    WV = permute(double(spk.waveforms), [1 3 2]);
    
elseif nargin == 3
    load(fn);
    T = double(spk.unwrapped_times);
    WV = permute(double(spk.waveforms), [1 3 2]);
    
    switch records_flag
        case 1 % implies that records_to_get is a timestamp list.
            keep = (T == records_to_get);
            T = T(keep);
            WV = WV(keep, :, :);
            
        case 2 %implies that records_to_get is a record number list
            keep = records_to_get;
            T = T(keep);
            WV = WV(keep, :, :);
            
        case 3 % implies that records_to_get is range of timestamps (a vector with 2 elements: a start and an end timestamp)
            keep = (T > records_to_get(1) & T <= records_to_get(2));
            T = T(keep);
            WV = WV(keep, :, :);
            
        case 4 % implies that records_to_get is a range of records (a vector with 2 elements: a start and an end record number)
            keep = records_to_get;
            T = T(keep(1):keep(2));
            WV = WV(keep(1):keep(2), :, :);
            
        case 5 % asks to return the count of spikes (records_to_get should be [] in this case)
            T = length(T);
    end
    
elseif nargin == 2
    % New "get" construction"
    if strcmp(a, 'get')
        switch (b)
            case 'ChannelValidity'
                T = [true true true true]; return;
            case 'ExpectedExtension'
                T = '.mat'; return;
            otherwise
                error('Unknown get condition.');
        end
    else
        error('2 argins requires "get" as the first argument.');
    end
end