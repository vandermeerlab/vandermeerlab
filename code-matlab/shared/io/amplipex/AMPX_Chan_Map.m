function [ map ] = AMPX_Chan_Map()
%AMPX_Chan_Map Quickly plot the layout of the probe with all the broken
%channels along with the best channels.
%

%%
run(FindFile('*Keys.m'))

layout88 = reshape(ExpKeys.Probe_layout, 8,8)';
layout64 = reshape(layout88, 1,64);
[~, all_best] = AMPX_BestChan(ExpKeys);
bests = [all_best.vm all_best.vl all_best.dm all_best.dl];
%%
figure(10)
maximize
for ichan = 1:64
    subtightplot(8,8, ichan)
    if isempty(intersect(layout64(ichan), bests)) ==0 && isempty(intersect(layout64(ichan), ExpKeys.BadChannels)) ==1; colour = [0 1 0];
    elseif isempty(intersect(layout64(ichan), ExpKeys.BadChannels)) ==0; colour = [1 0 0];
    else
        colour = [0 0 0];
    end
    text(.05, .05, num2str(layout64(ichan)), 'color', colour, 'FontSize', 36)
    xlim([ 0.025 .10]); ylim([0.025 .10])
    axis off
    
end

map = layout64;
