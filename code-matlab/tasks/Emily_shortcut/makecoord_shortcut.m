function [ coord ] = makecoord_shortcut(pos_x, pos_y, trajectory)

max_dist = 1;
newX = nan;
newY = nan;

for val = 1:length(trajectory)
    newX(end + 1) = trajectory(val, 1);
    newY(end + 1) = trajectory(val, 2);
end

coord = [newX(2:end); newY(2:end)];

% find the distance between each Coord point
diff_coord = diff(coord');
length_coord = sum((diff_coord.*diff_coord)').^0.5;

% Set up an array for the interpolated points
newPoints = [];
for iF = 1:length(length_coord)
    newPoints = [newPoints, coord(:, iF)];
    if length_coord(iF) > max_dist
        m = diff_coord(iF, 2)/diff_coord(iF, 1);
        b = (coord(2, iF) - m*(coord(1, iF)));
        nPoints = floor(length_coord(iF)/max_dist);  % number of points to add
        StepSize = -(coord(1, iF) - coord(1, iF+1))/nPoints;
        for iN = 1:nPoints - 1
            newY = m*(coord(1, iF) + StepSize*(iN)) + b;
            coord(1, iF);
            newPoints = [newPoints, [(coord(1, iF) + StepSize*(iN)) newY]'];
        end;
    end;
end
newPoints = [newPoints, coord(:,end)];
coord = newPoints;
end
