function [tf] = LoadExpKeys2
%LOADEXPKEYS2 Load ExpKeys into caller workspace
%   No need for any inputs; it checks the current directory automatically
%
%      tf = LoadExpKeys2
%        tf is 1 if the keys file was found and 0 if not
%
%        To suppress tf output, use [~] = LoadExpKeys2
%
% This function assumes that the keys file contains only one struct called
% 'ExpKeys'.
%
% A.Carey Feb, 2015

%%
[~,name,~] = fileparts(pwd); 

ugh = [name(1:4),'_',name(6:9),'_',name(11:12),'_',name(14:15),'_keys.m'];

fn = FindFiles(ugh);

if isempty(fn)
       tf = 0;     
else
    tf = 1;
    run(ugh)
    assignin('caller','ExpKeys',ExpKeys)
end

end

