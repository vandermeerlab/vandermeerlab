%%
in_fd{1} = 'D:\data\R120\R120-2017-05-20_pre-CNO';
in_fd{2} = 'D:\data\R120\R120-2017-05-20_post-CNO';

out_fd = 'D:\data\R120\R120-2017-05-20_pre-CNO';
out_fd_postfix = '_merged';

%% get filenames
for iFD = 1:length(in_fd)
    
   cd(in_fd{iFD});
   tt{iFD} = sort(FindFiles('*.ntt'));
   csc{iFD} = sort(FindFiles('*.ncs'));
    
end

%% match filenames
for iTT = 3:length(tt{1,1})
    
    fn1 = tt{1,1}{iTT};
    [~,fn1,~] = fileparts(fn1);
    
    for ifn2 = 1:length(tt{1,2})
        [~,fn2{ifn2},~] = fileparts(tt{1,2}{ifn2});
    end
    
    match_id = strmatch(fn1,fn2);
    
    if isempty(match_id)
        fprintf('No match found for %s...\n',fn1);
        continue;
    end
    
    % match found, do loading
    cd(in_fd{1});
    [Timestamps1, ScNumbers1, CellNumbers1, Features1, Samples1, Header1] = Nlx2MatSpike(tt{1,1}{iTT}, [1 1 1 1 1], 1, 1, [] );
    
    cd(in_fd{2});
    [Timestamps2, ScNumbers2, CellNumbers2, Features2, Samples2, Header2] = Nlx2MatSpike(tt{1,2}{match_id}, [1 1 1 1 1], 1, 1, [] );
    
    % merge
    Timestamps3 = cat(2,Timestamps1,Timestamps2);
    ScNumbers3 = cat(2,ScNumbers1,ScNumbers2);
    CellNumbers3 = cat(2,CellNumbers1,CellNumbers2);
    Features3 = cat(2,Features1,Features2);
    Samples3 = cat(3,Samples1,Samples2);
    
    % write
    out_fn = cat(2,fn1,out_fd_postfix,'.ntt');
    
    cd(out_fd);
    Mat2NlxSpike(out_fn, 0, 1, [], [1 1 1 1 1], Timestamps3, ScNumbers3, CellNumbers3, Features3, Samples3, Header1);
    
end