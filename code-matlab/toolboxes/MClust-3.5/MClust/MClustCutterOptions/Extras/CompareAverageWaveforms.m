function [redraw,rekey,undoable] = CompareAverageWaveforms()

% plots avergae waveforms of only shows

redraw = false; rekey = false; undoable = false;

global MClust_Clusters MClust_Hide

figure
iW = 1;
iH = 1;
toggle = 0;

showIt = ~MClust_Hide(1:(length(MClust_Clusters)+1));
num_clusters = sum(showIt);
iClust = 1;
while iW * iH < num_clusters
	if showIt(iClust)
		if toggle == 0
			iW = iW + 1;
			toggle = 1;
		else
			iH = iH + 1;
			toggle = 0;
		end
	end
	iClust = iClust+1;
end
plotCount = 1;
for iClust = 1:length(MClust_Clusters)
	if showIt(iClust)
		h = subplot(iH,iW,plotCount);
		ShowAverageWaveform(iClust,h);
		plotCount = plotCount+1;
	end
end

