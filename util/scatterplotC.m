function A=scatterplotC(x,y,C,varargin)
%  H=scatterplotC(x,y,C)
%  H=scatterplotC(x,[],C)
%
% Plots x and/or y with color proportional to C. 
%
%INPUTS
%  x,y -- input vector(s) to plot
%  C   -- color level to plot
%         (C can be an index 0 to N, where N is the number of colors)
%
%OUTPUTS
%  H   -- figure handle
%
%Parameters\Defaults
%
%  plotchar   = 'o'
%  NumColors  = 100;
%  Method     = 'Linear';
%  solid_face = 0;
% 
%  Note that all line graphic properties as in plot are supported

%
% JCJ 25 Oct 2002
% JCJ 01 May 2005 -- enabled varargin for use of all line graphic properties as in plot
%                 -- also included ability to make solid markers
%
% MvdM 14 Feb 2010  -- added fixed scale capability

plotchar   = 'o';
NumColors  = 100;
Method     = 'Linear';
LineWidth  = 1;
solid_face       = 0;

keepVarArgIn = [];
for iV = 1:2:length(varargin) % JCJ 01 May 2005 -- enabled varargin for use of all line graphic properties as in plot
	switch lower(varargin{iV})
		case 'solid_face'     %                 -- also included ability to make solid markers
			solid_face  = varargin{iV+1};
		case 'plotchar'
			plotchar  = varargin{iV+1};
		case 'numcolors'
            NumColors = varargin{iV+1};
		case 'method'
            Method    = varargin{iV+1};
        case 'scale'
            Scale = varargin{iV+1};
		otherwise
			keepVarArgIn=[keepVarArgIn iV iV+1];
	end
end
varargin = varargin(keepVarArgIn);


if exist('Scale','var')
    disp('Scale mode');
    extract_varargin;
    [garb,iC] = histc(C,linspace(Scale(1),Scale(2),NumColors));
elseif sum(round(C)==C)==length(C) % if color input is an index
    iC = C;
    NumColors=max(iC);
    extract_varargin;
else
    extract_varargin;
    if strncmpi(Method ,'Log',3)
        iC = floor(log(C)/log(max(C))*(NumColors-1))+1;
    elseif strncmpi(Method ,'Linear',3)
        iC = floor(C/max(C)*(NumColors-1))+1;
    else
        disp('Warning: Unrecognized Method. Using Linear. (scatterplotC)')
        iC = floor(C/max(C)*(NumColors-1))+1;
    end
end

if (length(iC)~=length(x))|((length(iC)~=length(y))&~isempty(y))
    disp('Warning: C has different length from x and/or y (scatterplotC)') 
end

caxis([0 NumColors]);
cmap = colormap(jet(NumColors));          % Make 'Hot' colormap for FR coding

H    = gcf;                               % Gets figure handle to plot on; if no figure, makes one
hold on

if ~isempty(y)
    for iColor = 1:NumColors
        FR2plot=find(iC==iColor); % Find array positions of spikes in firing rate range
        if ~isempty(FR2plot)
            if solid_face         %                 -- also included ability to make solid markers
                plot(x(FR2plot),y(FR2plot),plotchar,'Color', cmap(iColor,:),'MarkerFaceColor',cmap(iColor,:),varargin{:})           % Plot spike in that range with same color
            else
                plot(x(FR2plot),y(FR2plot),plotchar,'Color', cmap(iColor,:),varargin{:})           % Plot spike in that range with same color
            end
        end
    end
else
    t=1:length(x);
    for iColor = 1:NumColors
        FR2plot=find(iC==iColor); % Find array positions of spikes in firing rate range
        if ~isempty(FR2plot)
            if solid_face         %                 -- also included ability to make solid markers
                plot(t(FR2plot), x(FR2plot),plotchar,'Color', cmap(iColor,:),'MarkerFaceColor',cmap(iColor,:),varargin{:})           % Plot spike in that range with same color
            else
                plot(t(FR2plot), x(FR2plot),plotchar,'Color', cmap(iColor,:),varargin{:})           % Plot spike in that range with same color
            end
        end
    end
end

hold off

if nargout
    A=H;
end

