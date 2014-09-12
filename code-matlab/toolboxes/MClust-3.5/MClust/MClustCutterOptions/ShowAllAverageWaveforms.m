function [redraw,rekey,undoable] = ShowAllAverageWaveforms()

global MClust_Clusters 

redraw = false; rekey = false; undoable = false;
for iClust = 1:length(MClust_Clusters)
	ShowAverageWaveform(iClust);
    %if iClust > 1
	set(gcf, 'WindowStyle', 'docked');
    %end
end

% MvdM
desktop=com.mathworks.mde.desk.MLDesktop.getInstance;
desktop.setGroupDocked('Figures',false);
desktop.setGroupMaximized('Figures',true); 

