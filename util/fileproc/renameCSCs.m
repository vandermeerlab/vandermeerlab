function renameCSCs(fname,varargin)
% function renameCSCs(fname,varargin)
%
% WARNING: assumes normal config, i.e. CSC1-4 is TT1, 5-8 TT2, etc..
%
% MvdM 2012

ttsub_map = 'dabc';
extract_varargin;

flist = FindFiles('*.ncs');

for f = 1:length(flist)
	
    [~,n,e] = fileparts(flist{f});
    
    % get matching tetrode number
	cscno = sscanf(n,'CSC%d');
    ttid = ceil(cscno/4); 
    
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