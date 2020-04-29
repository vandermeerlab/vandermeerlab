%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                              ExpKeys                                %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ExpKeys.experimenter = {'aacarey'};
ExpKeys.species = 'Rat';
ExpKeys.behavior = 'MotivationalT';
ExpKeys.target = {'dCA1'};
ExpKeys.hemisphere = {'right'};

ExpKeys.day = 5;

ExpKeys.RestrictionType = 'food'; % 'food' % note, this is what was *withheld*, so rat is motivated to get this
ExpKeys.RestrictionDuration = []; % (hours) time from beginning of restriction to the beginning of experiment
ExpKeys.RatWeight = 581; % weight in grams, includes weight of array
ExpKeys.TimeWeighed = '11:30 am'; % time when the rat was weighed (began handling for experiment)


ExpKeys.Session = 'standard'; % 'pit', 'reversal’
ExpKeys.Layout = 'foodLeft'; % 'foodRight'
ExpKeys.Pedestal = 'L'; % did he start on the Left(L) or the Right(R) pedestal?
ExpKeys.pathlength = 334; % in cm; distance from start to finish (center arm to reward L or R). Note: start platform length was divided by two, because rat often landed on midpoint
ExpKeys.patharms = 369; % in cm; distance from reward to reward (L arm to R arm)
ExpKeys.realTrackDims = [185 167]; % x width and y width (according to camera axes)
ExpKeys.convFact = [2.6853 2.6169]; % x conversion factor and y conversion factor for loading position data in cm

ExpKeys.nPellets = 5;
ExpKeys.waterVolume = [];

ExpKeys.nTrials = 20;
ExpKeys.forcedTrials = [10 15 18]; % IDs of trials with block at choice point
ExpKeys.nonConsumptionTrials = []; %all consumed
ExpKeys.badTrials = 1; % rat attempted to reverse directions % IDs of trials where rat was somehow interfered with (e.g. unplanned block)

ExpKeys.TimeOnTrack = 5879; % with small buffer before "recording" segment begins
ExpKeys.TimeOffTrack = 8057; %with small buffer after "recording" segment ends
ExpKeys.prerecord = [4371.442104 5872.369604]; % timestamps%[1 3001856]; % [start stop] indices for prerecord, in terms of ntt tvec
ExpKeys.task = [5879.163104 8056.698604]; % timestamps %[3001857 7356928]; % [start stop] indices for task recording segment, in terms of ntt tvec
ExpKeys.postrecord = [8062.073104 8976.760604]; % timestamps %[7356929 9186304]; % [start stop] indices for postrecord, in terms of ntt tvec


ExpKeys.goodSWR = {'R050-2014-04-02-CSC07a.ncs' 'R050-2014-04-02-CSC08a.ncs' 'R050-2014-04-02-CSC06a.ncs'}; % list CSCs (at most one per TT) with good SWRs, first best
ExpKeys.goodTheta = {'R050-2014-04-02-HS1R2.ncs' 'R050-2014-04-02-HS3R2'}; % list CSCs (at most one per TT) with good theta, first best


