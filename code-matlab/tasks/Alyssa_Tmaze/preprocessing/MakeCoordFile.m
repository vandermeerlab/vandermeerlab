%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%      MAKE COORD FILE (DRAW TRAJECTORIES) for motivational shift     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script makes metadata.coord containing trajectories with
% coordinates in units of both pixels and centimeters. You need the version
% of MakeCoord that has varargins YDir and rot.
% A coord file contains the idealized trajectory that your rat would take
% along your track (rather than the wobbly one he takes in reality). Do not
% draw your coords to perfectly lie on the 90 degree bends in your track: 
% try to draw a good average of the rat's trajectory. So, if he cuts 
% corners, draw your coords with rounded edges.
%
% convFact is also created here, so add it to the corresponding ExpKeys script.
% Calculating the pixel-to-centimeter conversion factor requires the field 
% ExpKeys.realTrackDims. How to measure this properly is explained in the
% section "convert units to cm" down below.
%
% This script will probably work with other rats for tasks with TWO
% trajectory options. However, make sure you find out the value of
% rotation in [rotation,~] = view when you plot the raw position data.
% This will help you orient the figure so that LEFT and RIGHT are
% intuitively plotted. (Note: choice point might not be applicable for
% certain layouts.)
%
% aacarey Jan 2015
% edit March 2016

%% Set your current directory to the session you want to work with

clear

% load position data with units in pixels:
cfg = [];
pos = LoadPos(cfg);

% find out which rat we're working with
[~,sessionID,~] = fileparts(pwd);
rat = sessionID(1:4);

% the rotation of the position data when plotted in a figure
rotation = 0;
%[rotation,el] = view; % look at default values for rotation and el (they are 90,90),we want to change rotation to get the proper view

if isequal(rat,'R042'); rotation = 270; end % rotate the view for R042's data so the T is upright when plotted

%% This is what your position data looks like; is the maze orientation as you prefer?
% getd says "pull out the x data from pos"

% plot x pixels vs y pixels, reverse Y axis and/or change view to see proper maze orientation

figure; plot(getd(pos,'x'),getd(pos,'y')); title('Default, as exists in data'); xlabel('x data'); ylabel('y data');
view(rotation,90); set(gca,'YDir','reverse'); xlabel('x data'); ylabel('y data');

if isequal(rat,'R042')
    title('Figure rotated and flipped, preserving original axes values');
else
    % R044, R050, R064
    title('Figure flipped, values preserved');
end

%% draw trajectories for L and R trials (output is in pixel units)
% figure will open that will prompt you to click/draw out a trajectory.
% press enter to commit your trajectory.

coordL = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw left trajectory','YDir','reverse','rot',rotation); % CoordL is in units of pixels
coordR = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw right trajectory','YDir','reverse','rot',rotation); % CoordR is in units of pixels

%% click on choice point

% plot coordL and coordR so I can see where they sit

figure; set(gca,'YDir','reverse'); view(rotation,90); hold on;

plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4);
plot(coordL(1,:),coordL(2,:),'ob'); plot(coordR(1,:),coordR(2,:),'og'); title('Click choice point; press enter');
maximize

% get data values where user clicks
[x,y] = ginput; 

plot(x,y,'or','MarkerSize',10,'LineWidth',4); pause(1); close

%% convert units to cm 
% we want to use the exact same trajectories, but have them exist with
% different units
% You need to know the dimensions of your track, but these dimensions must
% be as if the track fits perfectly into a box that has the same orientation
% as the camera's field of view (the box has to be orthogonal to the field
% of view, or whatever)

% rhombus track ex: you don't want lengths a and b, you want x and y (I think, right?)

%    ..... x ..... 
%  _ _ _ _ _ _ _ _ _ 
% |                  |
% |        .         |  .
% |      .   .  a    |  .
% |    .       .     |  .
% |  .           .   |  y 
% |    .       .     |  .
% |      .   .  b    |  .
% |        .         |  .
% |                  |
% |_ _ _ _ _ _ _ _ _ |
%   field of view

% The T-maze is already approx lined up with the camera's field of view, so x = a and y = b
%     ......x.......
%  _ _ _ _ _ _ _ _ _ _ 
% |         a          |
% |    ____________    | .             
% |   |     |      |   | .
% |   |     |      |   | .
% |   |     |      |   | y
% |         | b        | .  
% |         |          | .
% |         |          | .
% |_ _ _ _ _ _ _ _ _ _ |

LoadExpKeys % to get real track dims

if ~isfield(ExpKeys,'convFact')
    convFact = PosCon(pos,ExpKeys.realTrackDims,'YDir','reverse');
else 
    convFact = ExpKeys.convFact;
end
  
%%%%%%% ****** ADD CONVFACT TO EXPKEYS FOR ALL SESSIONS  ******  %%%%%%%
%%

coordL_cm = coordL; % copy coordL under a new variable name, and apply some changes:
coordL_cm(1,:) = coordL_cm(1,:)./convFact(1); % apply x conversion
coordL_cm(2,:) = coordL_cm(2,:)./convFact(2); % apply y conversion

coordR_cm = coordR; % as above, for R instead
coordR_cm(1,:) = coordR_cm(1,:)./convFact(1); % apply x conversion
coordR_cm(2,:) = coordR_cm(2,:)./convFact(2); % apply y conversion


% convert choice point to cm
chp = [x; y];
chp_cm = [x/convFact(1); y/convFact(2)];

% NOTE when using coord_cm, you must also use LoadPos in cm, or your units
% are not the same!

% put it all in a struct for tighter packing in the base workspace (when loading variables later)
coord = struct('coordL',coordL,'coordL_cm',coordL_cm,'coordR',coordR,'coordR_cm',coordR_cm,'chp',chp,'chp_cm',chp_cm);

%% plot cm coord and choice point

figure; set(gca,'YDir','reverse'); view(rotation,90); hold on; 

title('Your coords in centimeters')
plot(coordL_cm(1,:),coordL_cm(2,:),'ob'); plot(coordR_cm(1,:),coordR_cm(2,:),'og'); plot(chp_cm(1),chp_cm(2),'or','MarkerSize',10,'LineWidth',4)

%% Save coord field in metadata

% WARNING: running this section overwrites existing coord field (if one exists already)!

% If metadata exists, it is loaded (LoadMetadata [at the time this was
% written] does not error if it does not find a metadata file). If metadata
% does not exist it is created here.
LoadMetadata
metadata.coord = coord;

% now save
save([sessionID,'-metadata.mat'],'metadata'); % this saves the specified variables under the given savename