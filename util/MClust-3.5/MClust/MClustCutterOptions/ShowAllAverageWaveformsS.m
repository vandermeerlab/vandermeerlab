function [redraw,rekey,undoable] = ShowAllAverageWaveformsS()

global MClust_Clusters 

redraw = false; rekey = false; undoable = false;
nClu = length(MClust_Clusters);

i = 0; j = 0;
while i*j < nClu
    if i > j
        j = j + 1;
    else
        i = i + 1;
    end
end

figure;

for iClust = 1:nClu
    subplot(i,j,iClust);
	ShowAverageWaveform(iClust,'PlotHere',1);
end