function [redraw,rekey,undoable] = ShowAllAverageWaveforms()

global MClust_Clusters 
figure
redraw = false; rekey = false; undoable = false;
iW = 1;
iH = 1;
toggle = 0;
num_clusters = length(MClust_Clusters);
while iW * iH < num_clusters
    if toggle == 0
        iW = iW + 1;
        toggle = 1;
    else
        iH = iH + 1;
        toggle = 0;
    end
end
for iClust = 1:length(MClust_Clusters)
    h = subplot(iH,iW,iClust);
    ShowAverageWaveform(iClust,h);
    %if iClust > 1
	%set(gcf, 'WindowStyle', 'docked');
    %end
end

% % MvdM
% desktop=com.mathworks.mde.desk.MLDesktop.getInstance;
% desktop.setGroupDocked('Figures',false); % heh
% desktop.setGroupMaximized('Figures',true); % heh

