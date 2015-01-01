function Coord = MakeCoord(x, y, varargin)
% function Coord = MakeCoord(x,y,varargin)
%
% Opens a figure and prompts the user to draw an idealized path. Returns an 
% array containing the x and y positions of points on the path. Commonly
% used for drawing paths or trajectories over position data.
% Important: if you have data from RR1 (UW) YDir should always be
% 'reverse', because the camera mirrors in the y axis.
%
% Each point is separated from the next by approximately 1 pixel (in practice it ranges
% from 1.0 to 1.1 pixels of separation between points)
%
% Output variable Coord should be saved in data folder for later use
%
% VARARGINS
%    titl: string; the figure title. Default: 'Select Linearized Path'
%    XDir / YDir: 'normal'(default) or 'reverse'; direction of increasing 
%                values on the specified axis
%    rot: number of degrees to rotate the image (clockwise!)     
%
% FOR NOTES ON CHANGING MAZE ORIENTATION 
% open MakeCoord and read the section titled "MAZE ORIENTATION EXAMPLES"
%
% original by NCST
% modified MvdM 08, 2014-06-24
% modified ACarey, 2014-12-31 (handles axes reversals and figure rotations)

%% MAZE ORIENTATION EXAMPLES

% MakeCoord by default plots your position data according to how the camera
% passes it to the aquisition system. But this might not be how the user imagines it. 
% If you recorded in RR1 at UW, the camera is mirrored in the y
% axis, so always pass in YDir as 'reverse' (UNLESS LoadPos is modified to change y automatically in the future). 
% You can play around with pos data using the plot() function to find out
% what varargins to pass into MakeCoord: 
% plot(pos_x,pos_y); set(gca,'YDir','reverse'); view(rot,90); % where rot is the
% degrees of rotation in the CLOCKWISE direction. (90 is the elevation;
% leave it as 90).

% REVERSING THE Y AXIS TO CORRECT MIRRORING:
% coord = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Omg upside down');
%  _ _ _ _ _ _ _ _ _ _ 
% |                    |  
% |         |          |
% |         |          |
% |         |          |
% |         |          |   
% |   L     |      R   |
% |   |     |      |   |
% |   |_____|______|   | 
% |                    |
% |_ _  _ _ _ _ _ _ _ _|
%      Camera view

% but in your head, the maze looks like this:
% coord = MakeCoord(getd(pos,'x'),getd(pos,'y'),'Ydir','reverse');
%  _ _ _ _ _ _ _ _ _ _ 
% |                    |
% |    ____________    |             
% |   |     |      |   |
% |   |     |      |   |
% |   L     |      R   |
% |         |          |   
% |         |          |
% |         |          |
% |_ _ _ _ _ _ _ _ _ _ |
%      In your head

% REVERSING THE Y-AXIS AND ROTATING BY 270 DEGREES
% Or maybe the camera sees this (note: L and R are wrong!):
% coord = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Omg rotated and mirrored');
%  _ _ _ _ _ _ _ _ _ _
% |                    |  
% |        R _ _ _     |
% |               |    |
% |               |    |
% |   ____________|    |   
% |               |    |
% |               |    |
% |        L _ _ _|    |
% |                    |
% |_ _  _ _ _ _ _ _ _ _|
%    Camera view

% but in your head, the maze looks like this:
% coord = MakeCoord(getd(pos,'y'),getd(pos,'x'),'YDir','reverse','rot',270);
%  _ _ _ _ _ _ _ _ _ _ 
% |                    |       
% |    ____________    |             
% |   |     |      |   |
% |   |     |      |   |
% |   L     |      R   |
% |         |          |   
% |         |          |
% |         |          |
% |_ _ _ _ _ _ _ _ _ _ |
%      In your head

%% 
MaxDist = 1; % Maximum separation between Coord points, used to linearly interpolate between user selected points.
newX = [];
newY = [];
titl = 'Select Linearized Path';
wraparound = 0;
YDir = 'normal'; % 'reverse' flips the y axis
XDir = 'normal'; % 'reverse' flips the x axis
rot = 0; % if user wants to change rotational orientation
extract_varargin;

if isempty(x)
    return
end

figure
plot(x,y,'.','Color',[0.7 0.7 0.7],'MarkerSize',4)
title(titl);
set(gca,'YDir',YDir,'XDir',XDir); view(rot,90); % view(az,el) zaxis is rot and elevation is 90.
xlabel('X data'); ylabel('Y data');
maximize;
hold on

if isempty(newX) || isempty(newY)
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
while any(diff(newX) == 0) || any(diff(newY) == 0)
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

% trying to get the function to plot how user wants to see it, but also
% output the coord correctly. So this hack seems to work (ACarey)
temp = Coord(2,:); 
Coord(2,:) = Coord(1,:); 
Coord(1,:) = temp; 
% end of hack

pause(2);
close 