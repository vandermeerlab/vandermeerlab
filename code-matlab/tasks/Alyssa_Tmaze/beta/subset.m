function Asub = subset(cfg_in,A,n)
%SUBSET Obtain uniformly distributed random subset for ts or iv data
%   Asub = SUBSET(cfg,A,n) randomly selects a subset of entries from A such
%   that Asub contains n entries. Entries are A.t for ts, and A.tstart and
%   A.tend for iv. usr data is properly accounted for in the subset process.
%
%    INPUTS
%         cfg: config struct with fields controlling function behavior
%           A: ts or iv struct
%           n: the number of entries you want in the subset.
%
%    OUTPUT
%        Asub: ts or iv struct containing a subset of the entries in A.
%
%    CONFIG OPTIONS
%       cfg.resetRNG = 1; If 1, temporarily resets the rng settings so that
%                different calls will produce the same subset. After running,
%                the original rng settings are reinstated. If 0, the current
%                rng state is used and repeated runs with the same cfg and
%                inputs will not produce the same output.
%       cfg.sort = 1; If 1, sort the entries in Asub so that they're in 
%                in temporally ascending order. If 0, don't sort.
%       cfg.verbose = 1; If 1, tell me how many went in and how many come
%                out; if 0, don't.
%
%          
% aacarey, Jan 2015. edit Dec 2015.

%% Set things up

cfg_def.resetRNG = 1; % reseed the rng or not
cfg_def.sort = 1; % sort output or not
cfg_def.verbose = 1; % talk to me or not

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

% identify which datatype we're working with
if isfield(A,'t') 
    datatype = 'ts';
    nEntries = length(A.t);
elseif isfield(A,'tstart')
    datatype = 'iv';
    nEntries = length(A.tstart);
else
    error('Input must be ts or iv datatype')
end

% make sure user isn't dumb
if n > nEntries
    error('Input already has fewer than n samples')
end

%% Match number of data entries

if cfg.resetRNG
    % get the current state of the rng
    s = rng;
    % seed the random number generator so that the output is predictable each
    % time the function is called on
    rng(1,'twister')
end

% shuffle randomly, then keep the first n samples
p = randperm(nEntries); % permutes a vector same length as A...uniform distribution
keep = p(1:n);

% return rng to original state
if cfg.resetRNG; rng(s); end

% sort if wanted
if cfg.sort; keep = sort(keep); end

%% Now select a subset and output stuff

cfg_temp.verbose = 0;

switch datatype
    case 'ts'
        Asub = SelectTS(cfg_temp,A,keep);
        if cfg.verbose; fprintf('%s: %d ts in, %d ts out\n',mfun,nEntries,length(Asub.t)); end
    case 'iv'
        Asub = SelectIV(cfg_temp,A,keep);
        if cfg.verbose; fprintf('%s: %d iv in, %d iv out\n',mfun,nEntries,length(Asub.tstart)); end
end

% keep a record of cfg history (note we are using A to bypass SelectTS or
% SelectIV from being added. If it's best to show these, then change it:
if isfield(Asub,'cfg')
    Asub.cfg.history.mfun = cat(1,A.cfg.history.mfun,mfilename);
    Asub.cfg.history.cfg = cat(1,A.cfg.history.cfg,{cfg});
    %Asub.cfg.history.mfun = cat(1,Asub.cfg.history.mfun,mfilename);
    %Asub.cfg.history.cfg = cat(1,Asub.cfg.history.cfg,{cfg});
else
    Asub.cfg.history.mfun = mfun;
    Asub.cfg.history.cfg = {cfg};
end

end


