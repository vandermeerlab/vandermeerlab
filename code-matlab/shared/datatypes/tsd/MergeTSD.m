function tsd_out = MergeTSD(cfg_in,TSD,varargin)
%MERGETSD Merge data from multiple TSDs into a single, new TSD
%   Corresponding data values from each TSD.data are combined according to
%   cfg.method. So, tsd1.data(1) is combined with tsd2.data(1) and so on.
%   If data values in one TSD are much larger than in another, you may
%   consider rescaling your TSDs beforehand so that one is not weighted
%   higher than another: see rescale() and rescmean().
%
% tsd_out = MERGETSD(cfg,tsd1,tsd2,...) 
%       Input the config and multiple [1x1 structs] of TSD datatype
%
% tsd_out = MERGETSD(cfg,TSD)
%       Input the config and a single [1xn struct] of TSD datatype
%
%   CONFIG OPTIONS
%       cfg.method = 'mean' (default)
%       'median' - the mean of the middle two data values in sorted order
%       'mean' - the average of the data values
%       'geometricmean' - the nth root of the mean of the data values
%         note: cfg.geometricmean is not compatible with TSDs that contain
%         negative values if the nth root is even 
%       'addition' - not implemented
%       'subtraction' - not implemented
%
%       cfg.verbose = 1 (default)
%           1 - Prints text to command window
%           0 - Does not print text to command window
%
%   INPUTS
%       cfg - config struct with field controlling function behaviour (see
%           CONFIG OPTIONS)
%       TSD - two possible formats
%             [1x1 struct] x n inputs (variable number of inputs)
%             [1xn struct] x 1 input
%
% Why use MergeTSD? 
%   1. To combine score vectors from different detectors
%   2. In some cases, you may want to merge your CSCs to get a general LFP. 
%
% aacarey Sept 2015
% aacarey edit Jan 2017 to handle single input [1xn struct]

cfg_def.method = 'mean';
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

switch nargin
    case 1
        error('Not enough input arguments')
    case 2
        if length(TSD) == 1
           error('Not enough TSD inputs') 
        end
        tvec = TSD(1).tvec;
        
    otherwise
        % make a 1xn struct
        tvec = TSD.tvec;
        Args = {TSD,varargin{:}};
        for iArg = 1:length(Args)
            if length(Args{iArg}.tvec) ~= length(tvec)
               error('All TSDs must contain the same number of samples') 
            end
            TSD(iArg) = Args{iArg};
        end        
end

for iTSD = 1:length(TSD)
    if ~ CheckTSD(TSD(iTSD))
        error('Inputs must be TSD.')
    end
end

% tell me what you are doing 
if cfg.verbose  
    fprintf('MergeTSD: outputting the %s of %d TSDs\n',cfg.method,size(TSD,2))  
end

% collect all the data fields into an array
TSD = struct2cell(TSD);
data = cell2mat(TSD(4,1,:));

% merge the data, as specified by user
switch cfg.method
    case 'median'
        data = median(data,3);
    case 'mean'
        data = mean(data,3);
    case 'geometricmean'
        data = nthroot(prod(data,3),numArgs);
    case 'sum'
        data = sum(data,3);
    otherwise
        error('Unrecognized cfg.method. Better check that spelling ^_^')
end

% form the output TSD

tsd_out = tsd(tvec,data);

% need to add label if csc???

tsd_out = History(tsd_out,mfun,cfg);

end

