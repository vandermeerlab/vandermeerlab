% loading
csc = LoadCSC('R002-2012-04-15-CSC06a.ncs');
S = LoadSpikes(FindFiles('*.t'));

% Load VT and NEV file
[VTTimeStamps, x, y] = Nlx2MatVT('VT1.nvt',[1 1 1 0 0 0],0,1,[]); VTTimeStamps = VTTimeStamps * 10^-6;
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV('Events.nev',[1 1 1 1 1],1,1,[]);

% TimeO
TimeOnTrack = 1860;
TimeOffTrack = 4310;

%% remove NaNs from VT -- NOTE should really fit a spline or kalman filter -- remove zeros
keeps = find(~isnan(x) & ~isnan(y) & VTTimeStamps > TimeOnTrack  & VTTimeStamps < TimeOffTrack);
VTTimeStamps = VTTimeStamps(keeps); x = x(keeps); y = y(keeps);

%% should restrict to track only


%% scatterfields
figure(1); clf
plot(x,y,'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
hold on; axis off;

for iC = 1:length(S)
   
    S{iC} = Restrict(S{iC},TimeOnTrack,TimeOffTrack);
    
    spk_x = interp1(VTTimeStamps,x,Data(S{iC}),'nearest');
    spk_y = interp1(VTTimeStamps,y,Data(S{iC}),'nearest');
    h = plot(spk_x,spk_y,'.','MarkerSize',5,'Color',[1 0 0]);
    pause;
    delete(h);
    
end

%% spectrogram
cscD = Data(csc); cscR=Range(csc); cscD=decimate(cscD,4); cscR = cscR(1:4:end);
Fs = 1./median(diff(cscR));

params.pad = 2; params.fpass = [0 250]; params.tapers = [4 0.5 2]; params.Fs = Fs;

cscD = locdetrend(cscD,params.Fs,[1 0.5]);
[Sg,T,F] = mtspecgramc(cscD,[0.5 0.1],params);

T = T + cscR(1);
imagesc(T,F,log10(Sg)'); axis xy

%% filter CSC

fband = [6 10];

tsdD = locdetrend(cscD,Fs,[1 0.5]);
[b,a] = cheby1(1,0.5,[fband(1)/(Fs/2) fband(2)/(Fs/2)]);
f = filtfilt(b,a,cscD);
cscF = tsd(cscR,f);

%plot(Range(cscF),Data(cscF))

%% extract phase (should really exclude low power first)
ph = hilbert(f);
ph = angle(ph);
cscP = tsd(cscR,ph);

%plot(Range(cscP),Data(cscP))
%hold on;

% plot histogram of spike phases
spk_ph = interp1(Range(cscP),Data(cscP),Data(S{3}));
hist(spk_ph,100);

%% phase scatterplot

figure(1); clf
plot(x,y,'.','MarkerSize',1,'Color',[0.5 0.5 0.5]);
hold on; axis off;

for iC = 1:length(S)
   
    spk_x = interp1(VTTimeStamps,x,Data(S{iC}),'nearest');
    spk_y = interp1(VTTimeStamps,y,Data(S{iC}),'nearest');
    spk_ph = interp1(Range(cscP),Data(cscP),Data(S{iC}));
    h = ScatterplotC(spk_x,spk_y,spk_ph,'solid_face',1);
    pause;
    delete(h);
    
end

%% get relevant times from events (Rob's code -- needs to be cleaned up & made into function)

indexes = strmatch('AD',EventStrings); %Erase 'AD packetloss...' lines
TTLs = TTLs(setdiff(1:length(TTLs),indexes));

zeroidx = strmatch('Null',EventStrings); % catch indexed zero pellet drops as recorded in events
zeroEvTS = EVTimeStamps(zeroidx); % get timestamps of zero pellet drops
if isempty(zeroEvTS); clear zeroEvTS; end;

EVTimeStamps = EVTimeStamps(setdiff(1:length(EVTimeStamps),indexes)); %Crop out AD PacketLoss from EVTimeStamps
%%
figure(2);
hold on;
set(gca,'YLim',[0 3]);
f1ts=zeros(size(TTLs)); f2ts=zeros(size(TTLs)); f3ts=zeros(size(TTLs));

nPellets = [1 3 5];

for k=1:length(TTLs)
	if TTLs(k)==1, 
		TTLlength=(EVTimeStamps(k+1)-EVTimeStamps(k));
		TTLlength=int32(TTLlength/700000);
		if TTLlength==nPellets(1),colour=[1 0 0]; f1ts(k)=(EVTimeStamps(k)); end;
		if TTLlength==nPellets(2),colour=[0 1 0]; f2ts(k)=(EVTimeStamps(k)); end;
		if TTLlength==nPellets(3),colour=[0 0 1]; f3ts(k)=(EVTimeStamps(k)); end;
		plot(EVTimeStamps(k),1,'.','Color',colour,'MarkerSize',7); % zeros will not show up in figure 2
	end;
	if TTLs(k)==2, 
		TTLlength=(EVTimeStamps(k+1)-EVTimeStamps(k));
		TTLlength=int32(TTLlength/700000);
		if TTLlength==nPellets(1),colour=[1 0 0]; f1ts(k)=(EVTimeStamps(k)); end;
		if TTLlength==nPellets(2),colour=[0 1 0]; f2ts(k)=(EVTimeStamps(k)); end;
		if TTLlength==nPellets(3),colour=[0 0 1]; f3ts(k)=(EVTimeStamps(k)); end;
		plot(EVTimeStamps(k),2,'.','Color',colour,'MarkerSize',7);
	end;
end

f1ts=f1ts(f1ts~=0); f1ts = f1ts*10^-6;
f2ts=f2ts(f2ts~=0); f2ts = f2ts*10^-6;
f3ts=f3ts(f3ts~=0); f3ts = f3ts*10^-6;
if exist('zeroEvTS','var'); f0ts=zeroEvTS(zeroEvTS~=0); f0ts=f0ts*10^-6; end

%% plot spikes around cue time
spikePETH(S{1},f1ts,'window',[-2 5],'dt',0.1);

%% Speed calculation instantaneous and averaged for reward drop types 

goodSamples = find(x ~= 0 & y ~= 0);
t = VTTimeStamps(goodSamples);
x = x(goodSamples);
y = y(goodSamples);

%t = t*10^-6; % convert to seconds

x = tsd(t',x'); y = tsd(t',y');
s=GetLinSpd(x,y);

if exist('f0ts','var'); f1ts=f0ts; end;

minspd = 1;
st=Data(s);
st=st(goodSamples);
[f1ts,f2ts,f3ts] = elimBadLaps (t,x,y,f1ts,f2ts,f3ts);

figure(3);
plot(1:length(Data(s)),Data(s));

[~,XOD1,SD1] = tsdPETH2(s,f1ts);
[~,XOD2,SD2] = tsdPETH2(s,f2ts);
[~,XOD3,SD3] = tsdPETH2(s,f3ts);

window = [-2 5];
dt = 0.25;
xw = window(1):dt:window(2);
t1 = f1ts(~isnan(f1ts));
t2 = f2ts(~isnan(f2ts));
t3 = f3ts(~isnan(f3ts));

nT1=length(t1);
nT2=length(t2);
nT3=length(t3);

	m1 = nanmean(XOD1); 
	m2 = nanmean(XOD2); 
	m3 = nanmean(XOD3); 
	se1 =  nanstd(XOD1)/sqrt(nT1+1);
	se2 =  nanstd(XOD2)/sqrt(nT2+1);
	se3 =  nanstd(XOD3)/sqrt(nT3+1);
figure(4);
	plot(xw,m1,'r',xw,m1+se1,'r:',xw,m1-se1,'r:',xw,m2,'g',xw,m2+se2,'g:',xw,m2-se2,'g:',xw,m3,'b',xw,m3+se3,'b:',xw,m3-se3,'b:');
	set(gca, 'XLim', window);
	ylabel('value (red 1 pellet, green 2 pellet, blue 3 pellet)')
	xlabel('peri-event (sec)');