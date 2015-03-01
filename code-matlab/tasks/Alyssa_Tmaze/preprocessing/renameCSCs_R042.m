function renameCSCs(fname,varargin)
% function renameCSCs(fname,varargin)
%
% WARNING: assumes normal config, i.e. CSC1-4 is TT1, 5-8 TT2, etc..
%
% MvdM 2012

tt_map = [9 9 9 9 ... % for each CSC, which is the corresponding tt
    10 10 10 10 ...
    11 11 11 11 ...
    12 12 12 12 ...
    13 13 13 13 ...
    14 14 14 14 ...
    15 15 15 15 ...
    16 16 16 16 ...
    1 1 1 1 ...
    2 2 2 2 ...
    3 3 3 3 ...
    4 4 4 4 ...
    5 5 5 5 ...
    6 6 6 6 ...
    7 7 7 7 ...
    8 8 8 8];

ttsub_map = 'dabc';
extract_varargin;

flist = FindFiles('*.ncs');

for f = 1:length(flist)
	
    [~,n,e] = fileparts(flist{f});
    
    % get matching tetrode number
	cscno = sscanf(n,'CSC%d');
    ttid = tt_map(cscno);
    
    % get subchannel
    ttsub = mod(cscno,4);
    ttsub = ttsub_map(ttsub+1);
    
    ttid = num2str(ttid);
    if length(ttid) == 1
        ttid = cat(2,'0',ttid);
    end

	f_out = cat(2,fname,'-CSC',ttid,ttsub,'.ncs');

    disp(sprintf('renaming %s%s to %s',n,e,f_out));
    
    command = cat(2,'! ren ',n,e,' ',f_out);
	evalc(command);
	
end