%% inspect raw data
fname = 'CSC20.ncs';
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

%%
plot(Timestamps,Samples(1,:),'k');
title('Click start and end times of HS detach events');
temp = ginput;
gap = temp(:,1);

t_start = cat(1,0,gap(2:2:end)); % t_start is for the intervals to keep
t_end = cat(1,gap(1:2:end),Inf);

%%
save('HS_detach_times.mat','t_start','t_end')
%%
fd = FindFiles('*.ntt');
fixTT12 = 0;

for iF = 1:length(fd)
    
    fname = fd{iF};
    
    disp(sprintf('Processing tt %s (%d/%d)...',fname,iF,length(fd)));
    
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(fname, [1 1 1 1 1], 1, 1, [] );
    
    if iF == 11 & fixTT12
        disp('***FIXING...');
        Samples(:,1,:) = Samples(:,1,:).*2; % if gain was set incorrectly
    end
    
    %% restrict
    keep_idx = [];
    for iR = 1:length(t_start)
        
        keep_idx = cat(2,keep_idx,find(Timestamps > t_start(iR) & Timestamps < t_end(iR)));
        
    end
    
    %% export
    fname_out = regexprep(fname,'\.','r\.');
    Mat2NlxSpike(fname_out, 0, 1, [], [1 1 1 1 1], Timestamps(keep_idx), ScNumbers(keep_idx), CellNumbers(keep_idx), Features(:,keep_idx), Samples(:,:,keep_idx), Header);
    
end