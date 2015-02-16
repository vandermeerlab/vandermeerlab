%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%               Loading and Linearizing Position Data                 %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Walkthrough/tutorial workflow

% FUNCTIONS
% LoadPos()
% getd()
% PosCon()
% MakeCoord()
% LinearizePos()

% SECTIONS:
% Loading and Plotting Position Data
% Getting Position Data Conversion Factors
% Making Coords
% Getting Choice Points Manually
% Linearizing Position Data

% BACKGROUND

% Position data is collected by an overhead camera and position tracking 
% algorithm as it views the bright, point source of light emitted by the 
% headstage LEDs. This position data is saved in the Neuralynx .nvt file.

% In the vandermeerlab codebase, position data is loaded using the function
% LoadPos(), and the x and y coordinate positions are accessed using the
% function getd().

% In analyses involving the use of trajectories, it is necessary to store
% the idealized path a rat would take along the track. These idealized
% paths, or coords, are saved for later use in scripts that order place 
% cells based on their field positions on the track.

% This workflow takes you through position data loading and linearizing
% using the example session R050-2014-04-02 recorded in RR1 at UW. R050's
% maze is T-shaped and has one choice point. The end of the left arm has a
% food reward and the end of the right arm has a water reward.

% A.Carey, Feb 2015

%% Current Directory

clear

% set your current directory to the session you want to work with:

directory = 'D:\data\R050\R050-2014-04-02'; 

cd(directory) 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%              Loading and Plotting Position Data                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% using an empty config loads the position data in the default units of pixels:
cfg = [];

% load the position data using LoadPos()
pos = LoadPos(cfg);

%% plot the position data

% position data is accessed using the function getd()
% if you want the x data points, x = getd(pos,'x');
figure; plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4); xlabel('x data'); ylabel('y data');
title('This maze is upside down!')

% In reality, R050's maze was oriented such that the top of the T points up. Notice
% here that the T-maze is actually upside down! This is because the camera
% mirrors the data in the Y axis.

% We can set certain axis properties so that the position data is
% plotted in the orientation we want to see it, but at the same time
% preserve the real data points and their relationships to one another:


%% Flip the Y axis upside down:

figure; plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4); xlabel('x data'); ylabel('y data');
set(gca,'YDir','reverse') % this flips the Y axis
title('This mazeis rightside up')


%% You can also rotate the maze:

figure; plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4); xlabel('x data'); ylabel('y data');
set(gca,'YDir','reverse') % this flips the Y axis
view(270,90); % the number 270 rotates the plot 270 degrees CLOCKWISE
title('Figure rotated AND flipped')

% Note that you can also perform this flip rotation by plotting (y,x), but
% this also transposes the data and works only for mazes that are bilaterally 
% symmetrical in the y axis (I think >_<). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     The above plotting examples are important for the next steps,   %%%
%%%     when making coords properly and getting the outer               %%%
%%%     boundaries of your track                                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%           Getting Position Data Conversion Factors                  %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BACKGROUND INFORMATION

% You might need to load position data in centimeters instead of in units
% of pixels.

% Do do this, you need the real dimensions of your track in cm but these 
% dimensions must be as if the track fits perfectly into a box that has the 
% same orientation as the camera's field of view (the box has to be orthogonal 
% to the field of view, or whatever word it is). See digrams below.

% Furthermore, the dimesions should be measured according to the trajectory
% the rat would have followed, ie: don't measure from one outer edge to the
% opposite outer edge, because the rat doesn't walk along the edge -- he
% walks along the center (unless he's an acrobat or weird).

% rhombus track ex: you don't want lengths a and b, you want x and y 

%    |.... x ....| 
%  _ _ _ _ _ _ _ _ _ 
% |                  |
% |        .         |  _
% |      .   .  a    |  .
% |    .       .     |  .
% |  .           .   |  y 
% |    .       .     |  .
% |      .   .  b    |  .
% |        .         |  _
% |                  |
% |_ _ _ _ _ _ _ _ _ |
%   field of view

% R050's T-maze is already approx lined up with the camera's field of view, so x = a and y = b
%     |.....x......|
%  _ _ _ _ _ _ _ _ _ _ 
% |         a          |
% |    ____________    | _             
% |   |     |      |   | .
% |   |     |      |   | .
% |   |     |      |   | y
% |         | b        | .  
% |         |          | .
% |         |          | _
% |_ _ _ _ _ _ _ _ _ _ |

% The dimensions should be saved in your ExpKeys as:
%           ExpKeys.realTrackDims = [xWidth yWidth];

%% Now get the conversion factors using the function PosCon()

% PosCon() will reorient your maze if you pass in certain varargins. 
% read the help documentation on PosCon for more intructions
% If you recorded in RR1 you should always flip the Y axis unless LoadPos() 
% is changed to load the Y data differently

realTrackDims = [185 167]; % x width and y width in centimeters

convFact = PosCon(pos,realTrackDims,'YDir','reverse');

% You should save convFact in ExpKeys
%          ExpKeys.convFact = [xConvFact yConvFact];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                    Making a Coord File                              %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open the function MakeCoord() and read its internal documentation for
% more information

% let's get the idealized trajectories; MakeCoord() takes varargins that can
% reorient your maze as you prefer to see it. If you recorded in RR1 you
% should always flip the Y axis unless LoadPos() is changed to do load the
% Y data differently

coordL = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw left trajectory, press enter when done','YDir','reverse'); % CoordL is in units of pixels
coordR = MakeCoord(getd(pos,'x'),getd(pos,'y'),'titl','Draw right trajectory, press enter when done','YDir','reverse'); % CoordR is in units of pixels

% these coords should be converted to units of cm using the convFact you already
% collected above:

coordL_cm = coordL; % copy coordL under a new variable name, and apply some changes:
coordL_cm(1,:) = coordL_cm(1,:)./convFact(1); % apply x conversion
coordL_cm(2,:) = coordL_cm(2,:)./convFact(2); % apply y conversion

coordR_cm = coordR; % as above, for R instead
coordR_cm(1,:) = coordR_cm(1,:)./convFact(1); % apply x conversion
coordR_cm(2,:) = coordR_cm(2,:)./convFact(2); % apply y conversion

% put it all in a struct for tighter packing in the base workspace (when loading variables later)
coord = struct('coordL',coordL,'coordL_cm',coordL_cm,'coordR',coordR,'coordR_cm',coordR_cm);

clear coordL coordL_cm coordR coordR_cm

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                   Getting Choice Points Manually                    %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If your maze has any choice points, it's probably a good idea to get the
% coordinates of those choice points. Here's how you can do this using a
% script:

figure; plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4); set(gca,'YDir','reverse'); hold on;
plot(coord.coordL(1,:),coord.coordL(2,:),'ob'); plot(coord.coordR(1,:),coord.coordR(2,:),'og'); title('Click choice point; press enter');
maximize

% get user input:
[x,y] = ginput;

plot(x,y,'or','MarkerSize',10,'LineWidth',4); pause(1); close

% convert choice point units

chp = [x; y];
chp_cm = [x/convFact(1); y/convFact(2)];

% add to coord
coord.chp = chp;
coord.chp_cm = chp_cm;

% coord can be saved as a metadata field, and metadata can be saved as a
% .mat file for later use

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                     Linearizing Position Data                       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here's a quick example of how to linearize the position data. You can
% think of linearizing as "pushing the position data onto the nearest point
% along the coord trajectory". 

% Normally you would have a separate script that defines the trial start
% and stop times; for the sake of this example, here are the trial
% intervals for all right trials that R050 did for the example session:

tstart = [5891.562065; 6717.689002; 7169.673936; 7634.938343; 7845.917153];
tend = [5983.521056; 6726.830367; 7187.524988; 7644.281508; 7856.827519];

trial_iv_R = iv(tstart,tend);

clear tstart tend

%% Load position data in units of centimeters

% do do this, you need to specify conFact as a config field
cfg.convFact = convFact;

pos = LoadPos(cfg);


%% Now, restrict the position data to the right trials only:

pos_R = restrict(pos,trial_iv_R);

% You can plot the output to see what this looks like:
figure; plot(getd(pos_R,'x'),getd(pos_R,'y')); set(gca,'YDir','reverse')

% note that the lines connecting the end of the right arm to the start of
% the track are not actually present in the data...they are just artifacts
% of plotting.

%% Linearize the position data 

% pass the coord in via the config
cfg.Coord = coord.coordR_cm;

% use LinearizePos()
pos_R_lin = LinearizePos(cfg,pos_R);

figure; plot(getd(pos_R_lin,'z'),getd(pos_R_lin,'z_dist')); xlabel('z'); ylabel('z dist');

% pos_R_lin.data(1), the z (or linearized position data), is the same
% regardless of whether you load position data in centimeters or pixels.
% pos_R_lin.data(2), the z dist (or the amount the real position was
% displaced when it was pushed onto z) is different depending on the units
% you choose when loading pos.
