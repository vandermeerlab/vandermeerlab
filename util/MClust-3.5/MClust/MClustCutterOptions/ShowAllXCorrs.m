function [redraw,rekey,undoable] = ShowAllXCorrs
% MvdM using ADR, cowen code

global MClust_TTData MClust_Clusters MClust_FeatureData MClust_Colors

redraw = false; rekey = false; undoable = false;
% get input parms
prompt={'bin size (msec):','Window width (msec):'};
def={'1','50'};
dlgTitle='ShowAllXCorrs';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);

binsize_msec = str2num(answer{1}); wwidth_msec = str2num(answer{2});
nC = length(MClust_Clusters);

fig_h = figure; set(fig_h,'Color',[0 0 0]);
spm = reshape(1:nC^2,[nC nC]); % to find subplot index for x,y plot

for c1 = 1:nC
	for c2 = c1:nC
			
		%ts1 = Range(ExtractCluster(MClust_TTData, FindInCluster(MClust_Clusters{c1}, MClust_FeatureData)),'ts');
		%ts2 = Range(ExtractCluster(MClust_TTData, FindInCluster(MClust_Clusters{c2}, MClust_FeatureData)),'ts');
        ts1 = Range(ExtractCluster(MClust_TTData, FindInCluster(MClust_Clusters{c1})),'ts');
		ts2 = Range(ExtractCluster(MClust_TTData, FindInCluster(MClust_Clusters{c2})),'ts');

		nbins = round(wwidth_msec/binsize_msec);
		[Y, xdim] = CrossCorr(ts1, ts2, binsize_msec, nbins);
		
		if (c1 == c2)
			col = MClust_Colors(c1+1,:);
			xlabel_text = num2str(c1);
			
			[maxval,maxind] = max(Y); % eliminate acorr = 1 at dt = 0
			Y(maxind) = 0;
		else
			col = [1 1 1];
			xlabel_text = [num2str(c1) '-' num2str(c2)];
		end
	
		sub_h = subplot(nC,nC,spm(c1,c2));
		b = bar(sub_h,xdim,Y);
		
		set(b,'FaceColor',col,'EdgeColor','none');
		axis tight; set(gca,'Color','none');
		%h = xlabel(xlabel_text); set(h,'FontSize',10,'Color',[1 1 1]);
		
	end
end


