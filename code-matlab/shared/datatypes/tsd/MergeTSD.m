function tsd_out = MergeTSD(cfg_in,varargin)
%MERGETSD Merge data from multiple TSDs into a single, new TSD
%   Corresponding data values from each TSD.data are combined according to
%   cfg.method. So, tsd1.data(1) is combined with tsd2.data(1) and so on.
%   If data values in one TSD are much larger than in another, you may
%   consider rescaling your TSDs beforehand so that one is not weighted
%   higher than another: see rescale() and rescmean().
%
% tsd_out = MergeTSD(cfg,tsd1,tsd2,...)
%
% cfg.method = 'mean'; 
%    'median' - the mean of the middle two data values in sorted order
%    'mean' - the average of the data values
%    'geometricmean' - the nth root of the mean of the data values
%         note: cfg.geometricmean is not compatible with TSDs that contain
%         negative values if the nth root is even 
%    'addition' - not implemented
%    'subtraction' - not implemented
%
% cfg.verbose = 1; % 1 - tell me what you are doing. 0 - don't 
%
% Why use MergeTSD? 
%   1. To combine score vectors from different detectors
%   2. In some cases, you may want to merge your CSCs to get a general LFP. 
%
% aacarey Sept 2015

cfg_def.method = 'mean';
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

numArgs = length(varargin);
if numArgs < 2
    error('Require at least two input TSDs')
end

% check that the inputs are TSDs
for iArg = 1:numArgs
    pass = CheckTSD(varargin{iArg});
    
    if ~ pass
       error([inputname(iArg+1),' must be a TSD object.'])
    end          
end

% check that the inputs all have the same number of samples 
checktvecs = nan(1,numArgs);
for iArg = 1:numArgs
     checktvecs(iArg) = length(varargin{iArg}.tvec);      
end
checktvecs = unique(checktvecs);
if length(checktvecs) > 1
    error('All TSDs must have the same number of samples.')
end

% tell me what you are doing 
if cfg.verbose
    names = [];
    for iArg = 1:numArgs
        if iArg < numArgs
            names = [names,' ',inputname(iArg+1)];           
        else
            names = [names,' and ', inputname(iArg+1)];
        end       
    end    
    disp(['MergeTSD: outputting the ',cfg.method,' of',names])   
end

% time to do the thing

% collect all the data fields into an array
allTSDs = nan(numArgs,length(varargin{iArg}.tvec));
for iArg = 1:numArgs
    allTSDs(iArg,:) = varargin{iArg}.data; 
end

% merge the data, as specified by user
switch cfg.method
    case 'median'
        data = median(allTSDs,1);
    case 'mean'
        data = mean(allTSDs,1);
    case 'geometricmean'
        data = nthroot(prod(allTSDs,1),numArgs);
    case 'sum'
        data = sum(allTSDs,1);
    otherwise
        error('Unrecognized cfg.method. Better check that spelling ^_^')
end

% form the output TSD

tsd_out = tsd(varargin{1}.tvec,data);

% need to add label if csc???

tsd_out = History(tsd_out,mfun,cfg);

end

