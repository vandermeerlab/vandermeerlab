%%% Remove VT mistakes and interpolate position data %%%

fn = FindFile('*.nvt');

% Load data
[Timestamps, X, Y, Angles, Targets, Points, Header] = Nlx2MatVT(fn, [1 1 1 1 1 1], 1, 1, [] );

% Identify VT mistakes 