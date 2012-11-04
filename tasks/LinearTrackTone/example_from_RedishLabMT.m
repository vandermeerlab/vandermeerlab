% MT single cell explorer. Run from directory above individual
% sessions.
%
% MvdM 07
%
% version 0

clear sd; close all;

%%%%%%%%%%%%%%%%
%%% SETTINGS %%%
%%%%%%%%%%%%%%%%

treemode = 0;

SET_WriteImage = 1;
SET_ReWriteImage = 1; % overwrite existing images

SET_RewardWindow = [0 5];
SET_RewardNiter = 1000;
SET_RewardWriteImage = 1;

SET_Recompute_AdaptiveFields = 0;
SET_Recompute_LinSpeedPETH = 0;

SET_FigBGCol = [1 1 1];

SET_CluFontSize = 10;
SET_CluMaxScale = 1000;
SET_CluMinScale = -500;

SET_AcorrLineWidth1 = 0.5; % unsmoothed
SET_AcorrLineWidth2 = 1; % smoothed
SET_AcorrLineWidth3 = 1; % baseline

SET_TextFontSize = 8;
SET_LabelsFontSize = 10;

SET_ISI2_Threshold = 2; % ISI cutoff for 2nd order plot

SET_PosSampleSize = 1;
SET_SpikeMarkerSize = 3;

SET_SamplingCheckBinRange = 200;

SET_RateMapNBins = 64;
SET_MinDwellTime = 0.2;
SET_SpeedSubSample = 1000; % subsample down to #spikes

%%%%%%%%%%%%%%%%%%%%
%%% LOADING DATA %%%
%%%%%%%%%%%%%%%%%%%%

if treemode

% get dir list to investigate
dirlist = dir; dirnames = {dirlist.name}; dirnames = dirnames([dirlist.isdir]); dirnames = dirnames(3:end);

% attempt to load
nSessions = length(dirnames);
nGoodSessions = 0;

else
	
	nSessions = 1;
	nGoodSessions = 0;
	dirnames{1} = '.';
	
end

for iS = 1:nSessions
%for iS = 1:1

    disp(sprintf('Inspecting session %s...',dirnames{iS}));
    
    % look at keys file to see what we have; look for MT unit data
    
	eval(['cd ' dirnames{iS}]);

	d = dir('*keys.m');
	
	if ~length(d)
		disp('No keys file found, Skipping...');
		cd ..
		continue;
	end
	
	[fd fn fe] = fileparts(FindFile('*keys.m'));
	disp(sprintf('Evaluating keys file %s...',fn));
	eval(fn);
	ExpKeys.fd = fd;
	
	if(strcmp(ExpKeys.Behavior,'MultipleT') && (length(ExpKeys.maze) == 4) && (ExpKeys.nPelletsPerFeeder == 2))
    
		disp(sprintf('Session %s is 4-T (2 pellets per feeder).',dirnames{iS}));
		
		% check if error file exists
		try
			%load(FindFile('*errors.mat'));
			%SessionData.errors = errors;
			FindFile('*errors.mat');
		catch
			disp('No errors file found. Eyeballing... -- hit key to continue --');
			EyeballErrors;
			pause;
		end
		
		if treemode
		cd ..
		end
		
		temp = MTinit(dirnames{iS});
		nCells = length(temp.S);
		
		if (length(temp.fc) ~= nCells)
			disp(sprintf('WARNING: nCells ~= length(fc) in session %s',dirnames{iS}));
		end
			
		if (nCells > 0)
			disp(sprintf('Session %s contains spiking data: added to data set.',dirnames{iS}));
			sd{nGoodSessions+1} = temp;
			goodDirNames{nGoodSessions+1} = dirnames{iS};
			nGoodSessions = nGoodSessions + 1;
		else
			disp(sprintf('No spiking data found in session %s, skipped.',dirnames{iS}));
		end
			
		% check if Coord file exists
		
		eval(['cd ' dirnames{iS}]);
		try
			load(FindFile('*Coord.mat'));
		catch
			disp('No Coord file found. Starting Raster_Lin_Coord...');
			Coord = Raster_Lin_Coord(temp.x,temp.y);
			save([fn(1:end-4) 'Coord'],'Coord');
		end
		if treemode
		cd ..
		end
		
		% check if IndMarker file exists
		
		eval(['cd ' dirnames{iS}]);
		try
			load(FindFile('*IndMarker.mat'));
		catch
			disp('No IndMarker file found. Starting Raster_Lin_IndMarker...');
			IndMarker = Raster_Lin_IndMarker(temp, Coord);
			save([fn(1:end-4) 'IndMarker'],'IndMarker');
		end
		if treemode
			cd ..
		end

		% check if speed file exists
		eval(['cd ' dirnames{iS}]);
		try
			load(FindFile('*_Speed.mat'));
		catch
			disp('No Speed file found. Generating...');
			spd = GetLinSpd(temp.x, temp.y);
			save([fn(1:end-4) 'Speed'],'spd');
		end
		if treemode
			cd ..
		end
		
		% check if adaptivefields exist
		eval(['cd ' dirnames{iS}]);
		try
			load(FindFile('*AdaptiveFields.mat'));
		catch
			disp('No AdaptiveFields file found. Generating...');
			[uac,af,afi,afout] = adaptive_field_script01(pwd);
			save([fn(1:end-4) 'AdaptiveFields'],'uac','af','afi','afout');
		end
		
		if SET_Recompute_AdaptiveFields
			disp('Regnerating AdaptiveFields...');
			[uac,af,afi,afout] = adaptive_field_script01(pwd);
			save([fn(1:end-4) 'AdaptiveFields'],'uac','af','afi','afout');
		end
		
		if treemode
			cd ..
		end
		
	else % if not a MultipleT session...
		if treemode
			cd ..
		end
	end
end


disp(sprintf('Starting analysis. Current directory: %s',pwd));

%%%%%%%%%%%%%%%%
%%% ANALYSIS %%%
%%%%%%%%%%%%%%%%

for iS = 1:nGoodSessions
%for iS = 1:1 % loop over sessions
    sess = sd{iS}; sessname = goodDirNames{iS};
	nCells = length(sess.S);
	
	% clear session vars
	clear CellData;
	
	eval(['cd ' sessname]);
	disp(sprintf('Session 1 (%s): just entered %s',sessname, pwd));
	
	
	[fd fn fe] = fileparts(FindFile('*keys.m'));
	disp(sprintf('Evaluating keys file %s...',fn));
	eval(fn);
	ExpKeys.fd = fd;
	
	% get error file
	try
		load(FindFile('*errors.mat'));
		%SessionData.errors = errors;
		sess.errors = errors;
	catch
		error('No errors file found.');
	end
	
	% Coord
	try
		load(FindFile('*Coord.mat'));
	catch
		error('No Coord file found.');
	end
	
	try
		load(FindFile('*IndMarker.mat'));
	catch
		error('No IndMarker file found.');
	end
	
	try
		load(FindFile('*_Speed.mat'));
	catch
		error('No Speed file found.');
	end
	
	try
		load(FindFile('*AdaptiveFields.mat'));
	catch
		error('No AdaptiveFields file found.');
	end
	
	% linearized pos
	zPos = GetLinPosition(data(sess.x), data(sess.y), Coord);
	zPos = tsd(Range(sess.x,'sec'), zPos);
	
	for iC = 25:25
		%for iC = 20:20 % loop over cells
	
		spks = sess.S{iC};
		cellfname = sess.fc{iC};
		
		if ~SET_ReWriteImage
			
			if (exist('CellProperties','dir'));
				cd CellProperties;

				fout = [sessname '_cell' num2str(iC) 'b.png']
					disp(sprintf('Image %s done already -- skipping.',fout));
					cd ..
					skipped = 1;
					
					if (exist(fout,'file'))
						try
							load(FindFile('*CellData.mat'));
						catch
							error('No CellData file found.');
						end
					
					continue;
				else
					skipped = 0;
				end

				cd ..
				
			end
		end
		
		%%%%%%%%%%%%%%
		%%% basics %%%
		%%%%%%%%%%%%%%
		
		nspks = length(Data(spks));
		tstart = StartTime(sess.x); tstop = EndTime(sess.x);
		ffreq = nspks/(tstop-tstart); % avg firing rate
		
		isis = diff(Data(spks));
		prop2 = sum(isis(isis > 2))/sum(isis); % Neil's PropISI>2
		ALL_prop2{iS}{iC} = prop2;
		[isih, binsUsed] = HistISI(spks);
		
		[Acorr1,binCenters] = AutoCorr(Data(spks)*10000,1,1000);
		Acorr1Conv = conv2(Acorr1,hamming(25),'same')./sum(hamming(25));
		pss = find(Acorr1Conv >= ffreq); pss = pss(1); % Neil's post-spike suppression
		ALL_pss{iS}{iC} = pss;
		
		[Acorr2,binCenters2] = AutoCorr(Data(spks)*10000,10,1000);
		Acorr2Conv = conv2(Acorr2,hamming(25),'same')./sum(hamming(25));
		
		%%%%%%%%%%%%%%%%%%%%%%
		%%% clustqual & wv %%%
		%%%%%%%%%%%%%%%%%%%%%%
		
		% find out filename base and cell number first
		%temp = regexp(cellfname,'(.*)_(\d+)\.','tokens'); % no idea why this works - .* surprisingly ungreedy?
		%clufname = [temp{1}{1} '-' temp{1}{2} '-ClusterQual'];
		%clufname = [regexprep(cellfname(1:end-2),'_','-') '-wv']; % for kendal
		clufname = [regexprep(cellfname,'_','-') '-wv'];
		[fp clufname fe] = fileparts(clufname);
		try
			eval(['load(''' clufname ''');']);
		catch
			try
				clufname = regexprep(clufname,'TT0','TT');
				eval(['load(''' clufname ''');']);
			catch
				try
					clufname = regexprep(clufname,'Sc','TT');
					eval(['load(''' clufname ''');']);
				catch
					error(sprintf('Unable to load wv file %s',clufname));
				end
			end
		end
		%wvfname = [temp{1}{1} '-' temp{1}{2} '-wv'];
		%wvfname = [regexprep(cellfname(1:end-2),'_','-') '-ClusterQual']; % for kendal
		wvfname = [regexprep(cellfname,'_','-') '-ClusterQual'];
		[fp wvfname fe] = fileparts(wvfname);
		try
			eval(['load(''' wvfname ''');']);
		catch
			try
				wvfname = regexprep(wvfname,'TT0','TT');
				eval(['load(''' wvfname ''');']);
			catch
				try
					wvfname = regexprep(wvfname,'Sc','TT');
					eval(['load(''' wvfname ''');']);
				catch
					error(sprintf('Unable to load CQ file %s',wvfname));
				end
			end
		end
		
		%%%%%%%%%%%%%%%%
		%%% FIGURE 1 %%%
		%%%%%%%%%%%%%%%%
		
		fig1h = figure(1); set(fig1h,'Color',SET_FigBGCol);
		sub1h = subplot('position',[0.05,0.55,0.4,0.4]);
		mWV = mWV / 10; sWV = sWV / 10;
		h = plot(xrange,mWV+sWV,'LineWidth',1,'Color',[0.5 0.5 0.5]); hold on;
		h = plot(xrange,mWV-sWV,'LineWidth',1,'Color',[0.5 0.5 0.5]);
		h = plot(xrange,mWV,'LineWidth',2);
		axis off;
		
		maxy = max(max(mWV+sWV))*1.2;
		axis([1 134 min(min(mWV-sWV)) maxy]);
		ampl = round(max(abs(max(mWV)-min(mWV))));
		
		h = text(0,maxy,sprintf('%s: %d spks (%.1f Hz)\nMaxAmpl: %d uV, ID %.1f, Lratio %f; PropISI>2 %.2f, PSS %d ms',cellfname,nspks,ffreq,ampl,CluSep.IsolationDist,CluSep.Lratio,prop2,pss), 'interpreter', 'none');
		set(h,'FontSize',SET_TextFontSize);
		
		%%%%%%%%%%%%
		%%% isih %%%
		%%%%%%%%%%%%
		sub2h = subplot('position',[0.55,0.7,0.4,0.25]);
		semilogx(binsUsed,isih);
		h = xlabel('isi (ms)'); set(h,'FontSize',SET_LabelsFontSize);
		set(gca,'FontSize',SET_LabelsFontSize);
		box off;
		
		%%%%%%%%%%%%%%
		%%% acorrs %%%
		%%%%%%%%%%%%%%
		
		sub2h = subplot('position',[0.05,0.35,0.4,0.1]);
		bar(binCenters(1:100),Acorr1(1:100));
		h = title('autocorrelations (scale in ms)'); set(h,'FontSize',SET_LabelsFontSize);
		set(gca,'FontSize',SET_LabelsFontSize);
		hold on; box off;
		%plot(binCenters(1:100),ffreq,'g','LineWidth',SET_AcorrLineWidth3);
		h = line([0 100],[ffreq ffreq],'Color',[0 0.5 0],'LineWidth',SET_AcorrLineWidth3);
		sub2h = subplot('position',[0.05,0.2,0.4,0.1]);
		plot(binCenters,Acorr1,'LineWidth',SET_AcorrLineWidth1);
		set(gca,'FontSize',SET_LabelsFontSize);
		hold on; box off;
		plot(binCenters,Acorr1Conv,'r','LineWidth',SET_AcorrLineWidth2);
		sub2h = subplot('position',[0.05,0.05,0.4,0.1]);
		plot(binCenters2,Acorr2,'LineWidth',SET_AcorrLineWidth1);
		set(gca,'FontSize',SET_LabelsFontSize);
		hold on; box off;
		%plot(binCenters2,Acorr2Conv(13:end-12),'r','LineWidth',SET_AcorrLineWidth2);
		
		%%%%%%%%%%%%%
		%%% frate %%%
		%%%%%%%%%%%%%
		sub2h = subplot('position',[0.55,0.45,0.4,0.15]);
		frate = FiringRate({spks},1);
		plot(Data(frate)); box off;
		h = xlabel('time (1s bins)'); set(h,'FontSize',SET_LabelsFontSize);
		h = ylabel('rate (Hz)'); set(h,'FontSize',SET_LabelsFontSize);
		set(gca,'FontSize',SET_LabelsFontSize);
		h = axis;
		axis([1 length(Data(frate)) h(3) h(4)]);
		
		%%%%%%%%%%%%%%%%%%%%%
		%%% 2nd order isi %%%
		%%%%%%%%%%%%%%%%%%%%%
		sub2h = subplot('position',[0.55,0.05,0.4,0.3]);
		isis_ = isis(isis < SET_ISI2_Threshold);
		%plot(isis_(1:end-1),isis_(2:end),'.');
		XData = isis_(1:end-1); YData = isis_(2:end);
		mkContours_simple(XData',YData',[0 SET_ISI2_Threshold],[0 SET_ISI2_Threshold],'Xnormalize',1);
		h = xlabel('isi(n) (s)'); set(h,'FontSize',SET_LabelsFontSize); h = ylabel('isi(n+1) (s)'); set(h,'FontSize',SET_LabelsFontSize);
		set(gca,'FontSize',SET_LabelsFontSize);
		axis([0 SET_ISI2_Threshold 0 SET_ISI2_Threshold]);
		
		%%%%%%%%%%%%%%%%%%%%%%%%
		%%% write image file %%%
		%%%%%%%%%%%%%%%%%%%%%%%%
		
		if SET_WriteImage
			
			if (~exist('CellProperties','dir'));
				mkdir('CellProperties');
			end
			cd CellProperties;
			fout = [sessname '_cell' num2str(iC)];
			eval(['print -f1 -dpng -r300 ' fout '.png']);
            eval(['print -f1 -depsc ' fout '.eps']);
			cd ..
			
		end

		close all;
		
		%%%%%%%%%%%%%%%%
		%%% FIGURE 2 %%%
		%%%%%%%%%%%%%%%%
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% position and spike position data %%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		[sfx,sfy] = ScatterFields({spks},sess.x,sess.y);
		
		fig2h = figure(2); set(fig2h,'Color',SET_FigBGCol);
		sub1h = subplot('position',[0.05,0.55,0.4,0.4]);
		plot(-Data(sess.y),-Data(sess.x),'.k','MarkerSize',SET_PosSampleSize);
		hold on; axis off;
		
		h = plot(-Data(sfy), -Data(sfx), 'm.','MarkerSize',SET_SpikeMarkerSize);
		
		%%% position sampling plot %%% 
		
		sub1h = subplot('position',[0.05,0.4,0.4,0.1]);
		bins = -SET_SamplingCheckBinRange:SET_SamplingCheckBinRange;
		y1 = hist(diff(Data(sess.x)),bins);
		y2 = hist(diff(Data(sess.y)),bins);
		semilogy(bins,mean([y1; y2]),'.'); box off;
		axis([bins(1)-0.5 bins(end)+0.5 0 max(mean([y1; y2]))]);
		set(gca,'FontSize',SET_LabelsFontSize);
		
		h = xlabel('dPixels'); set(h,'FontSize',SET_LabelsFontSize);
		h = ylabel('count'); set(h,'FontSize',SET_LabelsFontSize);
		h = title(cellfname,'interpreter','none'); set(h,'FontSize',SET_LabelsFontSize);
		
		%%% rate map %%%
		
		sub2h = subplot('position',[0.55,0.55,0.4,0.4]);
		[tc,occ] = TuningCurves(spks,sess.x,SET_RateMapNBins,sess.y,SET_RateMapNBins);
		
		occ(find(occ <= SET_MinDwellTime)) = NaN;
		pcolor(fliplr(flipud(tc{1}./occ)));
		axis off; shading flat; h = colorbar; set(h,'FontSize',SET_LabelsFontSize);
		
		%%% feeder PETHs %%%
		
		dt = 0.1;
		[outputS, outputIT] = spikePETH_simple(spks,sess.F1,'dt',dt);
		
		sub2h = subplot('position',[0.05,0.25,0.4,0.1]);
		if(length(outputS) > 1)
			m = histc(outputS, outputIT);
			bar(outputIT,m/dt/length(sess.F1)); shading flat;
            ax1 = gca;
        end
        
        % plot rate-bar graph
        axis([-2 5 0 35]); box off;
        set(gca,'FontSize',SET_LabelsFontSize,'LineWidth',1,'XTick',-2:5,'XTickLabel',-2:5,'YTick',0:10:30); box off;
        h = ylabel('rate (Hz)'); set(h,'FontSize',SET_LabelsFontSize);
        h = xlabel('time (s) from feeder 1 fire');  set(h,'FontSize',SET_LabelsFontSize);
        
        % plot speed
        ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','r','YColor','r');
		
        [P,X0D] = tsdPETH(spd,sess.F1,'dt',dt);
        hold on;
        plot(outputIT,data(P)/2.9,'Color','r','Parent',ax2);
        set(gca,'FontSize',SET_LabelsFontSize,'LineWidth',1,'XTick',-2:5,'XTickLabel',-2:5,'YTick',0:10:50); box off;
        axis([-2 5 0 50]);
        h = ylabel('speed (cm/s)'); set(h,'FontSize',SET_LabelsFontSize);
		
        [outputS, outputT] = spikePETH_simple(spks,sess.F2,'dt',dt);
		
		sub2h = subplot('position',[0.05,0.075,0.4,0.1]);
		if(length(outputS) > 1)
			m = histc(outputS, outputIT);
			bar(outputIT,m/dt/length(sess.F2)); shading flat;
            ax1 = gca;
        end
        
        % plot bar
        axis([-2 5 0 35]); box off;
        set(gca,'FontSize',SET_LabelsFontSize,'LineWidth',1,'XTick',-2:5,'XTickLabel',-2:5,'YTick',0:10:30); box off;
        h = ylabel('rate (Hz)'); set(h,'FontSize',SET_LabelsFontSize);
        h = xlabel('time (s) from feeder 2 fire');  set(h,'FontSize',SET_LabelsFontSize);
        
        % plot speed
        ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','bottom',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','r','YColor','r');
        
        [P,X0D] = tsdPETH(spd,sess.F2,'dt',dt);
        hold on;
        plot(outputIT,data(P)/2.9,'Color','r','Parent',ax2);
        set(gca,'FontSize',SET_LabelsFontSize,'LineWidth',1,'XTick',-2:5,'XTickLabel',-2:5,'YTick',0:10:50); box off;
        h = ylabel('speed (cm/s)'); set(h,'FontSize',SET_LabelsFontSize);
        axis([-2 5 0 50]);
		
		
		
		%%% speed PETH
		
		%[s,i] = LinSpdSTA(spks,spd);
		
		try
			load(FindFile('*LinSpeedPETH.mat'));
		catch
			disp('No LinSpeedPETH file found. Generating...');
			clear LinSpeedPETH;
			
			sTemp = Data(spks);
			
			if (length(sTemp) > SET_SpeedSubSample) % subsampling routine
				iTemp = randperm(length(sTemp));
				sTemp = sTemp(iTemp(1:SET_SpeedSubSample));
			end
			[linspeedPETH{iC},i] = tsdPETH(spd,sTemp);
			
		end
		
		sub2h = subplot('position',[0.55,0.075,0.35,0.35]);
		h = plot(Range(linspeedPETH{iC}),Data(linspeedPETH{iC}),'k','LineWidth',2);
		set(gca,'FontSize',SET_LabelsFontSize); axis tight; box off;
		h = xlabel('time (s) spike'); set(h,'FontSize',SET_LabelsFontSize);
		h = ylabel('linear speed (pix/s)'); set(h,'FontSize',SET_LabelsFontSize);
		
		if SET_WriteImage
			
			if (~exist('CellProperties','dir'));
				mkdir('CellProperties');
			end
			cd CellProperties;
			fout = [sessname '_cell' num2str(iC)];
			eval(['print -f2 -dpng -r300 ' fout 'b.png']);
            eval(['print -f2 -depsc ' fout 'b.eps']);
			cd ..
			
		end

		close all;
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% FANCY PLOTS FIGURE %%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        SET_DistanceRestrict = 1;
        SET_FilterWidth = 7; % error yes/no mode filter width (in samples)
        SET_dThr = 25; % distance to linpath for finding errors
        
        if SET_DistanceRestrict
            tvecE = Range(sess.x);
            xd = Data(sess.x); yd = Data(sess.x);
            tvecE = tvecE(~isnan(xd) & ~isnan(yd));

            [NN,d] = GetLinPosition(Data(sess.x),Data(sess.y),Coord);
            db = double(d > SET_dThr);
            db = medfilt1(db,SET_FilterWidth); % this works out to be a mode filter, because errors are {0 1}
            db = round(db);

            error_on = tvecE(find(diff([0; db]) == 1));
            error_off = tvecE(find(diff([0; db]) == -1));

            if length(error_on) == length(error_off)+1 % session ended in error?
                error_on = error_on(1:end-1);
            elseif length(error_on)+1 == length(error_off) % session started in error?
                error_off = error_off(1:end-1);
            elseif length(error_on) ~= length(error_off)
                error('Automatic error episode detection length mismatch.');
            end

            L0 = cat(2,sess.ExpKeys.TimeOnTrack,sess.F2);
            L1 = cat(2,sess.F2,sess.ExpKeys.TimeOffTrack);

            nonerror_on = cat(1,sess.ExpKeys.TimeOnTrack,error_off);
            nonerror_off = cat(1,error_on,sess.ExpKeys.TimeOffTrack);

            for iS = 1:length(sess.S)
                sessR.S{iS} = Restrict(sess.S{iS},nonerror_on,nonerror_off);
            end
            zPosR = Restrict(zPos,nonerror_on,nonerror_off);

        end
        sessR.F1 = sess.F1; sessR.F2 = sess.F2; sessR.errors = sess.errors;
        
		fig3h = figure(3); set(fig3h,'Color',SET_FigBGCol);
		sub1h = subplot('position',[0.075,0.75,0.8,0.2]); % linearized tuning curve
		RasterLinPlot(iC,zPosR,sessR,IndMarker,sub1h);
		
		sub2h = subplot('position',[0.075,0.075,0.4,0.6]); % lap raster
		RasterLinPlot(iC,zPosR,sessR,IndMarker,sub2h,'whattoplot',[0 0 0 0 1 0]);
		
		sub3h = subplot('position',[0.575,0.375,0.4,0.3]); % adaptive fields
		imagesc(squeeze(af(:,:,iC)));
		
		if SET_WriteImage
			
			if (~exist('CellProperties','dir'));
				mkdir('CellProperties');
			end
			cd CellProperties;
			fout = [sessname '_cell' num2str(iC)];
			eval(['print -f3 -dpng -r300 ' fout 'c.png']);
            eval(['print -f3 -depsc ' fout 'c.eps']);
			cd ..
			
		end
		
		
		
		%%% clean up, make celldata variable
		CellData{iC} = Parse_fn([cellfname '.t']);
		temp_tetrode = CellData{iC}.Tetrode; % weird MATLAB bug!
		CellData{iC}.depth = ExpKeys.TetrodeDepths(str2num(temp_tetrode));
		CellData{iC}.propisi = prop2;
		CellData{iC}.pss = pss;
		CellData{iC}.nspks = nspks;
		CellData{iC}.ampl = ampl;
		CellData{iC}.Lratio = CluSep.Lratio;
		CellData{iC}.IsolationDist = CluSep.IsolationDist;
		CellData{iC}.block = ExpKeys.block;
			
		%pause;
		close all;
		disp(sprintf('Cell %d done. NOTE: exiting, edit file to undo...',iC));
        return;
		
	end % loop over cells
	
	%try
	%	load(FindFile('*Reward.mat'));
	%catch
		disp('No Reward file found. Generating...');
		
		[RewardType,Perc1,Perc2] = rewardCell(sess,'window',SET_RewardWindow,'niter',SET_RewardNiter,'write_image',SET_RewardWriteImage);
		save([fn(1:end-4) 'Reward'],'RewardType','Perc1','Perc2','SET_RewardWindow','SET_RewardNiter');
						
	%end
		
	if exist('linspeedPETH','var')
		save([fn(1:end-4) 'LinSpeedPETH'],'linspeedPETH');
		%cd ..
	end
	
	if ~exist([fn(1:end-4) 'CellData'],'file')
			save([fn(1:end-4) 'CellData'],'CellData');
	end
	
	%if skipped
		cd ..
	%end
	
	disp(sprintf('Session %d done: now in dir %s',iS,pwd));
	
end % loop over sessions