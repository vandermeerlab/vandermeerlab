function [TSD,oldmin,oldmax] = RescaleTSD(TSD,newmin,newmax)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[TSD.data,oldmin,oldmax] = rescale(TSD.data,newmin,newmax);


end

