function renameTTs(fin,varargin)
% function renameTTs(varargin)
%
% example: renameTTs('R169-2009-06-26');
%
% MvdM 2007

TT1 = '01';
TT2 = '02';
TT3 = '03';
TT4 = '04';
TT5 = '05';
TT6 = '06';
TT7 = '07';
TT8 = '08';
TT9 = '09';
TT10 = '10';
TT11 = '11';
TT12 = '12';
TT13 = '13';
TT14 = '14';
TT15 = '15';
TT16 = '16';

fcount = 16;

extract_varargin;

for f = 1:fcount
	
    try
        fname = FindFile(cat(2,'TT',num2str(f),'.ntt'));
    catch
        disp(sprintf('No TT%d found, skipping',f));
        continue;
    end
    
	ttpart = eval(cat(2,'TT',num2str(f)));
	f_out = cat(2,fin,'-TT',ttpart,'.ntt');
	
    disp(sprintf('Renaming %s to %s',fname,f_out));
    
	command = cat(2,'! ren ',fname,' ',f_out);
	evalc(command);
	
end