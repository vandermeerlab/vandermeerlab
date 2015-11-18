%% MAKE COORD FILE (DRAW TRAJECTORIES) for T-maze motivational shift

% This script creates the coord.mat file containing trajectories with
% coordinates in units of both pixels and centimeters. You need the version
% of MakeCoord that has varargins YDir and rot.
% A coord file contains the idealized trajectory that your rat would take
% along your track (rather than the wobbly one he takes in reality)
% The coord information is useful later, in scripts that order place cells
% based on their positions on the track. 

% ACarey, Jan 1, 2015. Happy New Year...yaaaaaay (-_-)

%% Set your current directory to the session you want to work with

clear

% load position data with units in pixels:

cfg = [];
pos = LoadPos(cfg);

%% This is what your position data looks like; is the maze orientation as you prefer?
% getd says "pull out the x data from pos"

% plot x pixels vs y pixels, reverse Y axis and/or change view to see proper maze orientation

% When you find the desired rotation, you will pass the first variable of 
% view() into MakeCoord as the varargin 'rot' (R042 only)

% R050 (and R044?) T-maze

figure; plot(getd(pos,'x'),getd(pos,'y')); title('Default, as exists in data'); xlabel('x data'); ylabel('y data');
figure; plot(getd(pos,'x'),getd(pos,'y')); set(gca,'YDir','reverse'); title('Figure flipped, values preserved'); xlabel('x data'); ylabel('y data');

% R042 T-maze is rotated 90 deg compared to R050 and R044, so let's give it the same orientation in all plots

%figure; plot(getd(pos,'x'),getd(pos,'y')); title('Default, as exists in data'); xlabel('x data'); ylabel('y data');
%figure; plot(getd(pos,'x'),getd(pos,'y')); 
% %[az,el] = view; % look at default values for az and el (they are 90,90),we want to change az to get the proper view
%view(270,90); set(gca,'YDir','reverse'); title('Figure rotated and flipped, preserving original axes values'); xlabel('x data'); ylabel('y data');

%% draw trajectories for L and R trials (output is in pixel units)
% figure will open that will prompt you to click/draw out a trajectory.
% press enter to commit your trajectory.

% R050, R044
coordL = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw left trajectory','YDir','reverse'); % CoordL is in units of pixels
coordR = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw right trajectory','YDir','reverse'); % CoordR is in units of pixels


% R042
%coordL = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw left trajectory','YDir','reverse','rot',270); % CoordL is in units of pixels
%coordR = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw right trajectory','YDir','reverse','rot',270); % CoordR is in units of pixels

%% click on choice point

% plot coordL and coordR so I can see where they sit

figure; plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4); set(gca,'YDir','reverse'); hold on;
plot(coordL(1,:),coordL(2,:),'ob'); plot(coordR(1,:),coordR(2,:),'og'); title('Click choice point; press enter');
maximize

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

%% plot cm coord and cp

figure;set(gca,'YDir','reverse'); hold on; title('Your coords in centimeters')
plot(coordL_cm(1,:),coordL_cm(2,:),'ob'); plot(coordR_cm(1,:),coordR_cm(2,:),'og'); plot(chp_cm(1),chp_cm(2),'or','MarkerSize',10,'LineWidth',4)

%% Save coord field in metadata

% WARNING: running this section overwrites existing coord field (if one exists already)!

% first check if metadata exists yet
[~,name,~] = fileparts(pwd); % pwd is your current folder, we just want its namepart

loaded = LoadMetadata2;

if ~loaded % then it doesn't exist yet, so make a metadata struct
    metadata.coord = coord;
else % it does, so add a new field
    metadata.coord = coord;
end

% now save
savename = strcat(name,'-metadata.mat'); % use the folder's name, but concatenate it with '-metadata'
save(savename,'metadata'); % this saves the specified variables under the given [save]name