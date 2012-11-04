function [redraw, rekey, undoable]=ShowXCorrPairs

% ShowXCorrPairs()
%
% INPUTS
%   
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/GeneralizedCutterOptions folder

% JCJ Dec 2004 from ADR 2003

global MClust_Clusters

redraw   = false;
rekey    = false;
undoable = false;

prompt={'Clusters to Compare (space separated list):','bin size (msec):','Window width (msec):'};
def={'','1','500'};
dlgTitle='X Corr';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
if ~isempty(answer)
    iClust1=str2num(answer{1});
    nC1=length(iClust1)
    if nC1<1
        iClust1=1:length(MClust_Clusters);
        nC1=length(iClust1);
    end
    if nC1>1
        iClust2=iClust1(2:end);
        nC2=length(iClust2);
        [dim1, dim2]=rectdim(nC1*nC2/2);
    else
        iClust2=setxor(iClust1,1:length(MClust_Clusters));
        nC2=length(iClust2);
        [dim1, dim2]=rectdim(nC2);
    end
    figure
    cnt=0;
    for iC1=1:nC1
        for iC2=iC1:nC2
            if iClust1(iC1)~=iClust2(iC2) 
                cnt=cnt+1;
                subplot(dim1,dim2,cnt);
                MClustXcorr(iClust1(iC1), iClust2(iC2),  str2num(answer{2}),  str2num(answer{3}));
                a=axis;
                axis([-str2num(answer{3})/2 str2num(answer{3})/2, a(3:4)])
            end
        end
    end
end


function [dim1, dim2]=rectdim(N,Shape);
% [dim1, dim2]=rectdim(N);
% [dim1, dim2]=rectdim(N,Shape);
%
% Generates the number of rows and columns needed to display N elements in
% 2-D. If dim1 ~= dim2 then dim1 < dim2
%
% if Shape = 'square' gets dimensions as close to square but still finds
%              rectangular factors closest to a square if they exist.
%
% JCJ 2003-Oct-07

dim1=-inf;
dim2=+inf;

%convert Shape to a flag
if ~exist('Shape','var')
    Shape = 0;
else
    Shape = strcmpi(Shape,'square');
end



f=sort(factor(N));

while (length(f)<2) | (diff([dim1 dim2])>min([dim1 dim2]))
    
    f=sort(factor(N));
    
    
    if length(f)==2
        dim1=f(1);
        dim2=f(2);
    elseif length(f)>2
        P=perms(f);
        P=unique(P,'rows');
        temp=zeros(size(P,1),2);
        for iP=1:floor(length(f)/2)
            temp=[prod(P(:,1:iP),2) prod(P(:,iP+1:end),2)];
            [dG,iG]=min(diff(temp,[],2));
            if iP==1 || dG<diff(GoodPair)
                GoodPair = temp(iG,:);
            end
            
        end
        GoodPair=sort(GoodPair);
        dim1=GoodPair(1);
        dim2=GoodPair(2);
        
        
    end
    
    if Shape
        if (length(f)<2) | (diff([dim1 dim2])>min([dim1 dim2]))
            dim1=floor(sqrt(N));
            dim2=ceil(sqrt(N));
        end
        
        
        while (dim1*dim2)<N
            if dim1<dim2
                dim1=dim1+1;
            else
                dim2=dim2+1;
            end
        end
        return
    end
    
    N=N+1;
end

