function [t,x] = runlendec(p_in)
% function [t,x] = runlendec(p_in)
%
% converts a p-matrix describing posterior decodining probability [nTimeBins x nXbins], such as returned by
% Decode_Bayes_1S_1D.m, to arrays t and x corresponding to the original
% p-matrix ("inverse histogram" or run-length decoding).
%
% then, linear regression can be performed on the [t,x] data
%
% uses rude.m from the MATLAB Central FileExchange
%
% http://www.mathworks.com/matlabcentral/fileexchange/6436-rude--a-pedestrian-run-length-decoder-encoder
%
% MvdM 2015-02-07 initial version

f = 100; % would be nice to find some less arbitrary method for this
p_in = round(p_in.*f);

nTbins = size(p_in,1);
nXbins = size(p_in,2);

t = [];
x = [];

for iBin = 1:nTbins % could probably find some way to avoid this loop
    
   x_temp = rude(p_in(iBin,:),1:nXbins);
   t_temp = iBin.*ones(size(x_temp));
    
   t = cat(2,t,t_temp);
   x = cat(2,x,x_temp);
   
end
