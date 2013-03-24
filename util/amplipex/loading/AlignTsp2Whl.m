function AlignTsp2Whl(fbasename,colorvec)

% USAGE:
%     AlignTsp2Whl(fbasename,colorvec)
% 
% convert .tsp file to .whl file. Requires the .meta file
% 
% INPUT:
% 'fbasename': the file base names ('fbasename.tsp', etc.)
% 'colorvec': a 2 value vector [colFront colRear] which defines which
% color from [R G B] is the front and rear LEDs. Default: [1 3]
% 
% Adrien Peyrache adapted from Antal Berényi, 2012
% corrected by Antal Berenyi, 2012

% get timestamps and positions from tsp file
tspdata=load([fbasename '.tsp']);

%get sync delay in dat file
%load([fbasename '.concat.mat']); %delay in sample
%get length of dat file
infoFile = dir([fbasename '.dat']);
filelength=infoFile.bytes/(512*2);
clear infoFile
%get start and end timestamp of dat file
fid=fopen([fbasename '.meta']);
tline= fgetl(fid);
while ischar(tline)
    try
    if strcmp(tline(1:20),'TimeStamp of the end')
        tline=tline(59:end);
        EndTimestamp=sscanf(tline,'%d',1);
    end
    catch end
    try
    if strcmp(tline(1:22),'TimeStamp of the start')
        tline=tline(61:end);
        StartTimestamp=sscanf(tline,'%d',1);
    end
    catch
    end
    try
    if strcmp(tline(1:9),'Number of')
        tline=tline(31:end);
        ChanNum=sscanf(tline,'%d',1);
    end
    catch end
    try
    if strcmp(tline(1:9),'File size')
        tline=tline(21:end);
        DatSize=sscanf(tline,'%lu',1);
    end
    catch end
    
    tline= fgetl(fid);
end
fclose(fid);
    
%Calculate file length from dat file sample rat
DatLength=DatSize/(ChanNum*2*20); %Dat file size in ms
TspLength=EndTimestamp-StartTimestamp;

colorvec = sort([2*colorvec-1 2*colorvec])+1;
%remove lines from tspdata which has the same timestamp as the previous
%does
repeatingts=find(tspdata(1:end-1,1)==tspdata(2:end,1))+1;
for r=size(repeatingts):-1:1
    tspdata(repeatingts(r),:)=[];
end

%Clean and interpolating pos values BEFORE time interpolation
%tspclean = CleanWhlForR(tspdata(:,colorvec));
%tspclean(isnan(tspclean))=NaN;
%tspdata(:,colorvec) = tspclean;


%interpolate to 1 kHz - computer clock
t=tspdata(1,1):tspdata(end,1);
%interpTsp1kHz = zeros(length(t),length(colorvec)+1);
interpTsp1kHz = zeros(length(t),7);
interpTsp1kHz(:,1)=t;
warning off %Doen't like NaN but does well
%interpTsp1kHz(:,2:end)=interp1(tspdata(:,1),tspdata(:,colorvec),interpTsp1kHz(:,1));
interpTsp1kHz(:,2:end)=interp1(tspdata(:,1),tspdata(:,2:7),interpTsp1kHz(:,1));
warning on

%Find segments which are -1s in tspdata
for i=1:3%length(colorvec)
    minusonesegments(:,1)=tspdata(find(tspdata((2:end-1),i*2)==-1),1);
    minusonesegments(:,2)=tspdata(find(tspdata((2:end-1),i*2)==-1)+2,1);
    for j=1:size(minusonesegments,1)
        interpTsp1kHz(minusonesegments(j,1)-interpTsp1kHz(1,1)+1:minusonesegments(j,2)-interpTsp1kHz(1,1),i*2:i*2+1)=-1;
    end
    clear minusonesegments
end

%align the beginning
if (StartTimestamp<tspdata(1,1)) 
    i=1:tspdata(1,1)-StartTimestamp;
    tspnew(i,1)=StartTimestamp+i-1;
    %tspnew(i,2:length(colorvec)+1)=NaN;
    tspnew(i,2:7)=-1;
    tspnew=[tspnew; interpTsp1kHz];
else
    tspnew=interpTsp1kHz(StartTimestamp-tspdata(1,1)+1:end,:);
end
clear i

%adjust the end
% if (DatLength<size(tspnew,1)) 
%     tspnew=tspnew(1:DatLength,:);
% else
%     t = tspnew(end,1)+1:tspnew(end,1)+DatLength-size(tspnew,1);
%     %tspnew=[tspnew; [t' repmat(NaN,[DatLength-size(tspnew,1) length(colorvec)])]];
%     tspnew=[tspnew; [t' repmat(NaN,[DatLength-size(tspnew,1) 6])]];
% end
% clear t interpTsp1kHz

if (EndTimestamp<tspnew(end,1)) 
    tspnew=tspnew(1:find(tspnew(:,1)==EndTimestamp),:);
else
    t = tspnew(end,1)+1:EndTimestamp;
    %tspnew=[tspnew; [t' repmat(NaN,[DatLength-size(tspnew,1) length(colorvec)])]];
    tspnew=[tspnew; [t' repmat(-1,[size(t,2) 6])]];
end
clear t interpTsp1kHz

t = tspnew(1,1):size(tspnew,1)/DatLength:tspnew(end,1);
%interpTsp = zeros(length(t),length(colorvec)+1);
interpTsp = zeros(length(t),7);
interpTsp(:,1)=tspnew(1,1):size(tspnew,1)/DatLength:tspnew(end,1);
warning off
interpTsp(:,2:end)=interp1(tspnew(:,1),tspnew(:,2:end),interpTsp(:,1));
warning on



if 0 %Do we want a whl1k file? Personnally, no...
    fid=fopen([file '.whl1k'],'w');
    for i=1:size(interpTsp,1)
        fprintf (fid,'%i %i\n',interpTsp(i,2:3));
    end
    disp(['Samples in .dat file per channel: ' int2str(DatLength)]);
    disp(['Lines in .whl1k file:' int2str(size(interpTsp,1))]);
    fclose(fid);
end
%resample to 39.0625 Hz
%whldata=resample(interpTsp(:,2:end),390625,10000000);
whldata=interpTsp(floor(1:25.6:size(interpTsp,1)),2:end);

if 1
    figure(1),clf
    scatter(whldata(:,1),whldata(:,2),'.')
    %keyboard
end

dlmwrite([fbasename '.whlall'],whldata,'delimiter','\t');
dlmwrite([fbasename '.whl'],whldata(:,colorvec-1),'delimiter','\t');

end
