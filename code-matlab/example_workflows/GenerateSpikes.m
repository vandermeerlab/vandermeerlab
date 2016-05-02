%% Example script for using genInhomogenousPoisson
% genInhomogenousPoisson uses a deletion method for generating poisson spikes so it's quite fast but
% the way it's written was specified for the replay generative model and sort of hacky since we were
% going back and forth between position based and time based algorithms. (see pfmodel2 line 222)
% 
% Anyways, here I'll show a simple example of how to use it (and if need be alter it) so that it can
% be used in other situations.
%
% youkitan 2016-05-01 first
%
%
%% TL;DR EXAMPLE

% specify your rate function. For example, a sine wave:
tvec = 0:.01:2; %some time vector
ratefunc = @(t) 5*sin(10*t); %arbitrary function

% check that the function works
time_rate = ratefunc(tvec);
figure; subplot(2,1,1)
plot(tvec,time_rate)

% initiate input variables
max_rate = max(time_rate); %get it from the rate function;
total_time = tvec(end)-tvec(1); %just get it from your tvec

% generate spikes
spiketimes = genInhomogeneousPoisson(max_rate,total_time,ratefunc);
S = ts; S.t = spiketimes; %convert to ts struct for easy plotting

% check that it make sense
subplot(2,1,2); PlotSpikeRaster2([],S); 
axis([tvec(1),tvec(end),0,2]) %align axis with rate function



%% IN DEPTH EXAMPLE (as used in pfmodel2)
% total_time needs to be a duration but it can be extracted from an interval...
t1 = 0.5;
t2 = 2.5;
total_time = t2-t1;

% ratefunc is any function that takes a 1xN time vector as an input and returns a 1xN rate vector.
% Originally it was designed so that we had a function which created rates based on a position
% vector such as a gaussian function...
ratefunc1 = @(x, a, b, c) a.*exp(-(x-b).^2/(2*(c.^2)));

% we can simulate a place field as such:
tvec = t1:.01:t2; % 
pvec = 5*tvec; % position moves at a constant speed of 5
center = 8; %center position of gaussian
std = 1; %width of gaussian
max_rate = 10; %this will set the hight of the rate function as you'll see
pos_rate = ratefunc1(pvec,max_rate,center,std);

% you can see here that the rate function has a gaussian relationship with position
figure; subplot(3,1,1)
plotyy(tvec,pvec,tvec,pos_rate);

% then to use the spike generator we revert it to a time based function. We can our defined rate
% function for position as follows...
ratefunc2 = @(t) ratefunc1(interp1(tvec,pvec,t),max_rate,center,std);

% or if you have a given position vector (e.g., position tsd) and a tuning curve, you can
% interpolate the rate at a given time
ratefunc2 = @(t) interp1(pvec,pos_rate,interp1(tvec,pvec,t));


% you can see here that it returns the same thing as the "position" version...
time_rate = ratefunc2(tvec);
subplot(3,1,2); plot(tvec,time_rate);

% this means that you can actually use empirical place field distributions and convert them into a
% time-based rate function. Now you can use the spike generator with this time-rate function
spiketimes = genInhomogeneousPoisson(max_rate,total_time,ratefunc2);
S = ts; S.t = spiketimes; %convert to ts struct for easy plotting

% viola!
subplot(3,1,3); h2 = PlotSpikeRaster2([],S); 
axis([0.5,2.5,0,2]) %align axis with rate function

%% BONUS

% you can also simulate large batches of spiking with the nTrials parameter

% if one is interested, the gamma-poisson distribution (commented out in the function) generates a
% similar distribution of spike but with greater variance (depending on the beta parameter)


