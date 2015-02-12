function Asub = predsubset(A,n)
%PREDSUBSET Predictable random subset (an oxymoron!)
%   Asub = predsubset(A,B) randomly selects a subset of samples from A such
%   that Asub contains n samples. The samples are selected using a rng seed
%   that creates a predictable output. A can be S (cell array of spiketrains) or 
%   an iv data type (specifically trial intervals).  
%
%   If A is S (cell array of spiketrains), predsubset chooses a subset of 
%   spiketrains for Asub; If A is iv (trials), then predsubset chooses a 
%   subset of trials for Asub.
%          
% ACarey, Jan 2015. Written for T-maze project. Testing phase only!

%% Check input type
%A = S_pc_R; 

if isfield(A,'t') 
    datatype = 1;
elseif isfield(A,'tstart')
    datatype = 2;
else
    error('Input must be S or iv')
end


%% Match number of samples

% pull out the relevant part of the input struct
if datatype == 1
    A_orig = A;
    A = A.t; 
else
    A_orig = A; % keep original for tend and cfg in switch
    A = A.tstart; 
end

% make sure user isn't dumb
check = length(A);
if n > check
    error('Input already has fewer than n samples')
end

% get the current state of the rng
s = rng;
% seed the random number generator so that the output is predictable each
% time the function is called on
rng(1,'twister')

% shuffle A pseudorandomly, then keep the first n samples
p = randperm(length(A)); % permutes a vector same length as A...uniform distribution
A_shuffled = A(p);
Asubset = A_shuffled(1:n);

%% Now return the correct datatype struct

switch datatype
    case 1
        Asub.t = Asubset;
        % also need to shuffle and select subset for rest of fields (ugly)
        Asub.label = A_orig.label(p);
        Asub.label = Asub.label(1:n);
        Asub.cfg = A_orig.cfg;
        Asub.usr.rating = A_orig.usr.rating(p);
        Asub.usr.rating = Asub.usr.rating(1:n);
    case 2
        Asub.tstart = Asubset;
        % also need to shuffle and select subset for .tend
        tend_shuffled = A_orig.tend(p);
        Asub.tend = tend_shuffled(1:n);
        Asub.cfg = A_orig.cfg;
end

% return rng to original state

rng(s)

end


