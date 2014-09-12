function Coord = MakeCoord(x, y, varargin)
% function Coord = MakeCoord(x,y,varargin)
%
% Returns an array containing the x and y positions of points on an
% idealized path, input by the user
%
% Each point is separated from the next by approximately 1 pixel (in practice it ranges
% from 1.0 to 1.1 pixels of separation between points)
%
% Output variable Coord should be saved in data folder for later use
%
% original by NCST
% modified MvdM 08, 2014-06-24

MaxDist = 1; % Maximum separation between Coord points, used to linearly interpolate between user selected points.
newX = [];
newY = [];
titl = 'Select Linearized Path';
wraparound = 0;
extract_varargin;

if isempty(x)
    return
end

figure
plot(y,x,'.','Color',[0.7 0.7 0.7],'MarkerSize',4)
title(titl);
maximize;
hold on

if isempty(newX) | isempty(newY)
	[newY,newX] = ginput;
	newX = abs(newX)';
	newY = abs(newY)';
end

% wrap around the path
if wraparound
    newX(end + 1) = newX(1);
    newY(end + 1) = newY(1);
end

% add 1 pixel of jitter to any consequtive points which have the same x or y values
while any(diff(newX) == 0) | any(diff(newY) == 0)
    newX(find(diff(newX) == 0)) = newX(find(diff(newX) == 0)) + 1;
    newY(find(diff(newY) == 0)) = newY(find(diff(newY) == 0)) + 1;
end

Coord = [newX; newY];

% find the distance between each Coord point
dCoord = diff(Coord');
lCoord = sum((dCoord.*dCoord)').^0.5;

% Set up an array for the interpolated points
newPoints = [];
for iF = 1:length(lCoord)
    newPoints = [newPoints, Coord(:,iF)];
    if lCoord(iF) > MaxDist
        m = dCoord(iF,2)/dCoord(iF,1);
        b = (Coord(2,iF) - m*(Coord(1,iF)));
        nPoints = floor(lCoord(iF)/MaxDist);  % number of points to add
        StepSize = -(Coord(1,iF) - Coord(1,iF + 1))/nPoints;
        for iN = 1:nPoints - 1
            newY = m*(Coord(1,iF) + StepSize*(iN)) + b;
            Coord(1,iF);
            newPoints = [newPoints, [(Coord(1,iF) + StepSize*(iN)) newY]'];
        end;
    end;
end
newPoints = [newPoints, Coord(:,end)];
Coord = newPoints;

% Show the idealized path
plot(Coord(2,:),Coord(1,:),'og');
plot(Coord(2,1),Coord(1,1),'*b');
pause(2);
close 