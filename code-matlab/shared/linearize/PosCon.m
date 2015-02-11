function convFact = PosCon(pos,realTrackDims,varargin)
%POSCON Get position data conversion factors
%   convFact = PosCon(pos,realTrackDims,varargin)
%
% Opens a figure and prompts the user to define outer X width and Y width
% track boundaries. Using two clicks, first define the outer boundaries of
% the track in the X dimension. When you have completed two clicks the
% width will plot (no need to hit Enter), and you will be asked to do the 
% same for the Y dimension. The figure closes automatically.
% 
% Outputs are the dimension-specific conversion factors needed to convert
% position data from pixels to centimeters in LoadPos(). Make sure these
% are placed in ExpKeys and input to LoadPos() via the config.
%     ExpKeys.convFact = convFact;
%     cfg_pos.convFact = ExpKeys.convFact;
% 
% Important: if you have data from RR1 (UW) YDir should always be
% 'reverse', because the camera mirrors in the y axis. (See MakeCoord()
% function's internal documentation for more information).
%
% Inputs:
%        pos - position data; default output from LoadPos() when convFact
%        field not specified. 
%        realTrackDims - [realXwidth realYwidth] in centimers; ideally the
%        dimensions are based on the trajectory the rat would follow along
%        the track, rather the track edges.
%
% Varargins
%        XDir = 'normal'; % 'reverse' flips the x axis 
%        YDir = 'normal'; % 'reverse' flips the y axis
%        rot = 0; % if user wants to change rotational orientation
%
% Outputs
%       convFact - [xConvFact yConvFact];
%          -> xConvFact - conversion factor for x data points in pos
%          -> yConvFact - conversion factor for y data points in pos
%
% A.Carey, Feb 2015.

%%

YDir = 'normal'; % 'reverse' flips the y axis
%YDir = 'reverse';
XDir = 'normal'; % 'reverse' flips the x axis
rot = 0; % if user wants to change rotational orientation
extract_varargin;

if isempty(getd(pos,'x'))
    return
end

figure; plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4)
title('Define X boundary with two mouse clicks');
set(gca,'YDir',YDir,'XDir',XDir); view(rot,90); % view(az,el) zaxis is rot and elevation is 90.
xlabel('X data'); ylabel('Y data');
maximize;
hold on

[x,yX] = ginput(2);

plot([x(1) x(2)],[yX(1) yX(2)],'LineWidth',4,'Color','g');

%pause(1);

title('Define Y boundary with two mouse clicks');

[xY,y] = ginput(2);

plot([xY(1) xY(2)],[y(1) y(2)],'LineWidth',4,'Color','b');

pause(1);
close

x = sort(x); % we don't want negative conversion factors, right?
y = sort(y);

xdiff_pixel = x(2) - x(1);
xConvFact = xdiff_pixel / realTrackDims(1); % conversion factor for x data
  
ydiff_pixel = y(2) - y(1);
yConvFact = ydiff_pixel / realTrackDims(2); % conversion factor for y data

convFact = [xConvFact yConvFact];

end

